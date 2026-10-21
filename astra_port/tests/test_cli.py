#!/usr/bin/env python3
"""Run the installed-style CLI without a scientific data directory."""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--program',type=Path,required=True)
    p.add_argument('--examples',type=Path,required=True)
    p.add_argument('--netcdf',action='store_true')
    p.add_argument('--reference',type=Path)
    a=p.parse_args()
    program=a.program.resolve()
    try:
        import numpy as np
        import netCDF4
    except ImportError:
        np=netCDF4=None
        if a.reference: raise
    compared=0
    worst=0.

    def same(actual,expected,scale=0):
        nonlocal compared,worst
        actual,expected=np.asarray(actual),np.asarray(expected)
        assert actual.shape==expected.shape,(actual.shape,expected.shape)
        assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
        norm=np.maximum(np.maximum(np.abs(actual),np.abs(expected)),scale)
        error=np.divide(np.abs(actual-expected),norm,out=np.zeros_like(norm),where=norm!=0)/np.finfo(float).eps
        compared+=actual.size
        worst=max(worst,float(error.max(initial=0)))
        assert np.all(error<=64),f'{error.max()} epsilon-scaled error'

    with tempfile.TemporaryDirectory(prefix='turbotuv-cli-') as work:
        work=Path(work)

        def run(*args,ok=True):
            result=subprocess.run([str(program),*map(str,args)],cwd=work,capture_output=True,text=True)
            assert (result.returncode==0)==ok,(args,result.returncode,result.stderr)
            return result

        assert 'TurboTUV' in run('--version').stdout
        assert 'Usage:' in run('--help').stdout
        for args in [(),('--unknown',),('--output',),('absent.json',),
                     (a.examples/'tuv_5_4.json','--zenith','nan'),
                     (a.examples/'tuv_5_4.json','--zenith','181')]:
            run(*args,ok=False)
        if not a.netcdf:
            run(a.examples/'tuv_5_4.json','--netcdf',work/'disabled.nc',ok=False)
            assert not (work/'disabled.nc').exists()
        for example in ('tuv_5_4','ts1_tsmlt'):
            json_path=work/'results.json'
            run(a.examples/(example+'.json'),'--output',json_path)
            actual=json.loads(json_path.read_text())
            yaml_result=json.loads(run(a.examples/(example+'.yml')).stdout)
            assert actual==yaml_result,'JSON/YAML scientific results differ'
            assert actual['format']=='TurboTUV-1'
            nz,nw=len(actual['altitude_km']),len(actual['wavelength_nm'])
            assert all(len(s['direct_flux'])==nw and len(s['direct_flux'][0])==nz for s in actual['steps'])
            if a.reference:
                for kind,filename in [('photolysis','photolysis_rate_constants.nc'),('heating','heating_rates.nc'),('dose','dose_rates.nc')]:
                    labels=actual[kind+'_labels']
                    if not labels: continue
                    with netCDF4.Dataset(a.reference/('example_'+example)/filename) as ref:
                        def read(name):
                            var=ref[name]
                            # The original append writer places time last.
                            return np.moveaxis(var[:],var.dimensions.index('time'),0) if 'time' in var.dimensions else var[:]
                        same(actual['altitude_km'],ref['altitude'][:])
                        same(actual['wavelength_nm'],ref['wavelength'][:])
                        same([s['time_hours'] for s in actual['steps']],ref['time'][:])
                        same([s['solar_zenith_angle'] for s in actual['steps']],ref['solar zenith angle'][:],90)
                        same([s['radiation_scale'] for s in actual['steps']],ref['Earth-Sun distance'][:])
                        same([s['temperature_K'] for s in actual['steps']],read('temperature'))
                        for i,label in enumerate(labels):
                            same([s[kind][i] for s in actual['steps']],read(label.replace('/','_')))
                        for name,key in [('direct radiation','direct_flux'),('downward radiation','down_flux'),('upward radiation','up_flux')]:
                            if name in ref.variables: same([s[key] for s in actual['steps']],read(name),1)
            if a.netcdf:
                nc=work/'results.nc'
                run(a.examples/(example+'.json'),'--netcdf',nc)
                assert nc.read_bytes()[:8]==b'\x89HDF\r\n\x1a\n'
                if netCDF4 is not None:
                    with netCDF4.Dataset(nc) as ds:
                        same(actual['altitude_km'],ds['altitude'][:])
                        same(actual['wavelength_nm'],ds['wavelength'][:])
                        for name,key in [('direct radiation','direct_flux'),('downward radiation','down_flux'),('upward radiation','up_flux'),
                                         ('direct irradiance','direct_irradiance'),('downward irradiance','down_irradiance'),('upward irradiance','up_irradiance'),
                                         ('actinic flux','actinic_flux'),('spectral irradiance','spectral_irradiance')]:
                            same([s[key] for s in actual['steps']],ds[name][:])
                        for kind in ('photolysis','heating','dose'):
                            for i,label in enumerate(actual[kind+'_labels']):
                                same([s[kind][i] for s in actual['steps']],ds.groups[kind][label.replace('/','_')][:])
            zenith=json.loads(run(a.examples/(example+'.json'),'--zenith','30','--scale','1').stdout)
            assert len(zenith['steps'])==1
            assert all(row[-1]==1 for row in zenith['steps'][0]['direct_flux'])
            dark=json.loads(run(a.examples/(example+'.json'),'--zenith','180').stdout)
            assert all(value==0 for row in dark['steps'][0]['photolysis'] for value in row)
        # Published files keep the replaced file's mode, or get the umask default.
        import os,stat
        mask=os.umask(0); os.umask(mask)
        fresh=work/'fresh.json'
        run(a.examples/'tuv_5_4.json','--zenith','30','--output',fresh)
        assert stat.S_IMODE(fresh.stat().st_mode)==0o666&~mask,oct(fresh.stat().st_mode)
        shared=work/'shared.json'
        shared.write_text('{}'); shared.chmod(0o640)
        run(a.examples/'tuv_5_4.json','--zenith','30','--output',shared)
        assert stat.S_IMODE(shared.stat().st_mode)==0o640,oct(shared.stat().st_mode)
        # A failed run must not leave a complete JSON document on stdout.
        hostile=work/'albedo.json'
        config=json.loads((a.examples/'tuv_5_4.json').read_text())
        for profile in config['profiles']:
            if profile['name']=='surface albedo': profile['uniform value']=2
        hostile.write_text(json.dumps(config))
        failed=run(hostile,'--zenith','30',ok=False)
        try:
            json.loads(failed.stdout); raise AssertionError('failed run printed complete JSON')
        except json.JSONDecodeError: pass
        sentinel=work/'preserved.json'
        sentinel.write_text('preserve existing file')
        bad=work/'bad.json'
        bad.write_text('{grids: [')
        run(bad,'--output',sentinel,ok=False)
        assert sentinel.read_text()=='preserve existing file'
        assert not list(work.glob('*.tmp.*')),'temporary output leaked'
        assert not (work/'data').exists(),'test unexpectedly supplied external scientific files'
    print(f'CLI: JSON/YAML, isolated runtime, output integrity, {compared} numeric comparisons; worst {worst:.3g} eps')

if __name__=='__main__': main()
