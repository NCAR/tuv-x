#!/usr/bin/env python3
"""Extract complete spectral constructor cases from the unmodified TUV-x inputs.

Unit-test configuration type labels sometimes differ from the constructor called
by their Fortran test. Only those labels are normalized to factory names here.
Scientific options and datasets are unchanged. No runtime source files needed.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def generate(source):
    cases = []
    # Grid configurations carry the name and units TUV-x requires.
    def grid(path): return {'name': 'wavelength', 'type': 'from csv file', 'units': 'nm', 'file path': path}
    combined = grid('data/grids/wavelength/combined.grid')
    simple = grid('test/data/grid.photo.simple.csv')
    la = grid('test/data/grid.photo.la_sr.csv')
    types = {'hcfc': 'HCFC+hv->Products', 'clono2': 'ClONO2', 'ch2o': 'CH2O',
             'o3_tint': 'O3', 'rayliegh': 'air', 'n2o-n2_o1d': 'N2O+hv->N2+O(1D)',
             'ch3ono2_ch3o_no2': 'CH3ONO2+hv->CH3O+NO2',
             'nitroxy': 'CH3COCH2(ONO2)+hv->CH3COCH2(O.)+NO2'}
    sources = (sorted((source/'test/data/cross_sections').glob('*.json')) +
               sorted((source/'test/data/quantum_yields').glob('*.json')) +
               sorted((source/'test/data/spectral_weights').glob('*.json')) +
               [source/'test/data/heating_rates.json', source/'test/data/photorates.test.config.json'] +
               sorted((source/'examples').glob('*.json')) +
               [source/'data/dose_rates.json', source/'data/photolysis_rate_constants.json'])
    def walk(node, grid, origin, type_override=None):
        if isinstance(node, list):
            for child in node: walk(child, grid, origin, type_override)
        elif isinstance(node, dict):
            grid = next((g for g in node.get('grids', []) if g.get('name') == 'wavelength'), grid)
            for key, child in node.items():
                role = {'cross section': 0, 'cross sections': 0, 'quantum yield': 1,
                        'quantum yields': 1, 'weights': 2}.get(key)
                children = child if isinstance(child, list) else [child]
                if role is not None:
                    for item in children:
                        if isinstance(item, dict) and 'type' in item:
                            item = dict(item)
                            if role == 0 and type_override: item['type'] = type_override
                            cases.append((role, origin, grid, item))
                else: walk(child, grid, origin, type_override)
    provenance = []
    for path in sources:
        raw = path.read_text()
        if path.parent.name == 'data' and path.name in ('dose_rates.json','photolysis_rate_constants.json'):
            raw = '{' + raw + '}'
        raw = raw.replace("#PHOTO_RATES", (source/"data/photolysis_rate_constants.json").read_text())
        data = json.loads(raw)
        origin = path.relative_to(source).as_posix()
        grid, override = combined, None
        if 'test/data/cross_sections/' in origin:
            label = path.name.split('.')[1]
            grid = la if label == 'base' else simple
            override = types.get(label)
        elif origin == 'test/data/quantum_yields/base.config.json': grid = la
        walk(data, grid, origin, override)
        provenance.append({'path': origin, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    # Additional reference cases exercise constructor interactions absent from
    # the upstream scenarios, while using the same compiled scientific tables.
    narrow = {'name': 'wavelength', 'type': 'equal interval', 'units': 'nm', 'begins at': 100,
              'ends at': 400, 'cell delta': 1.25}
    for method in ('conserving', 'fractional source', 'fractional target', 'linear'):
        for fold in (False, True) if 'fractional' in method else (False,):
            file = {'file path': 'test/data/cross_sections/cross_section.base.nc',
                    'interpolator': {'type': method, 'fold in': fold},
                    'zero below': 101.25, 'zero above': 320,
                    'lower extrapolation': {'type': 'constant', 'value': 0.125}}
            cases.append((0, 'additional/interpolation', narrow, {'type': 'base', 'netcdf files': [file]}))
            cases.append((0, 'additional/merge', narrow, {'type': 'base', 'merge data': True,
                'netcdf files': [file, {'file path': 'test/data/cross_sections/cross_section.base.nc',
                                      'zero below': 102.5, 'zero above': 250}]}))
    for lower in ({'type':'boundary'}, {'type':'constant','value':0.123}):
        cases.append((1, 'additional/yield-upper-boundary', narrow, {'type':'base',
            'netcdf files':['test/data/quantum_yields/base.nc'],
            'lower extrapolation':lower, 'upper extrapolation':{'type':'boundary'}}))
    cases.append((0, 'additional/existing-point-data', la, {'type':'base',
        'netcdf files':[{'file path':'test/data/cross_sections/cross_section.base.nc'}],
        'data':{'default value':42, 'point values':[{'wavelength':102.5,'value':92.3}]}}))
    for basename in ('taylor','burkholder'):
        path = source/f'test/data/cross_sections/util/{basename}.config.json'
        param = json.loads(path.read_text())
        cases.append((0, f'additional/{basename}', combined, {'type':'temperature based','parameterization':param}))
        provenance.append({'path':path.relative_to(source).as_posix(), 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    param_grid = {'type':'equal interval','name':'wavelength','units':'nm',
                  'begins at':200,'ends at':360,'cell delta':1}
    for log in ('base 10','natural'):
        for inverse in (False,True):
            for with_file in (False,True):
                param = {'AA':[-20,0.001], 'BB':[0.0001,-0.000001], 'lp':[0,1],
                    'base temperature':298,'base wavelength':250,'logarithm':log,
                    'invert temperature offset':inverse, 'minimum wavelength':205,
                    'maximum wavelength':340, 'temperature ranges':[
                        {'maximum':220,'fixed value':220}, {'minimum':210,'maximum':350},
                        {'minimum':350,'fixed value':350}]}
                config = {'type':'temperature based','parameterization':param,
                          'parameterization wavelength grid':param_grid}
                if with_file: config['netcdf file']='data/cross_sections/H2O2_1.nc'
                cases.append((0,'additional/exponential-merge',combined,config))
    harwood = {'type':'HARWOOD','aa':[-20,-21],'bb':[1,2],'base temperature':100,
               'base wavelength':300,'logarithm':'natural','invert temperature offset':True}
    cases.append((0,'additional/harwood',combined,{'type':'temperature based',
        'parameterization':harwood,'parameterization wavelength grid':{
            'type':'from config file','units':'nm','name':'wavelength','values':[300,350,400]}}))
    for coeff in ([0.2], [-2,0.02,-0.00002]):
        cases.append((1,'additional/taylor-override',combined,{'type':'Taylor series',
            'constant value':0.1,'coefficients':coeff,'minimum wavelength':200,
            'maximum wavelength':330, 'override bands':[{'band':'range',
                'minimum wavelength':300,'maximum wavelength':400,'value':0.33}]}))
    blocks = ['/* Generated by tools/generate_spectrum_cases.py; do not edit. */\n',
              'static const struct spectrum_case { int role; const char *origin, *grid, *config; } spectrum_cases[] = {\n']
    for role, origin, grid, config in cases:
        fields = [str(role), json.dumps(origin), json.dumps(json.dumps(grid, separators=(',',':'))),
                  json.dumps(json.dumps(config, separators=(',',':')))]
        blocks.append('    {' + ', '.join(fields) + '},\n')
    blocks.append('};\n')
    return ''.join(blocks), json.dumps(provenance, indent=2)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT.parent/'TUV-x')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    content = generate(args.source)
    for name, text in zip(('spectrum_cases.inc','spectrum_cases_manifest.json'), content):
        path = ROOT/'tests/generated'/name
        if args.check:
            if path.read_text() != text: raise SystemExit(f'Stale file: {path}')
        else: path.write_text(text)

if __name__ == '__main__': main()
