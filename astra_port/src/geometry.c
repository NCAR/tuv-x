#include "turbotuv.h"
#include <math.h>
#include <stdint.h>

tt_status tt_spherical_geometry(size_t n, const double *z, double zenith, int *nid, double *dsdh) {
    if (n < 2 || n > INT32_MAX || !z || !nid || !dsdh || !isfinite(zenith) || zenith < 0.0 ||
        zenith > 180.0)
        return TT_INVALID;
    for (size_t i = 0; i < n; ++i)
        if (!isfinite(z[i]) || (i && z[i] <= z[i - 1]))
            return TT_INVALID;
    const double radius = 6.371e3, pi = 3.1415926535898;
    double re = radius + z[0];
    if (re <= 0.0)
        return TT_RANGE;
    double sinrad = sin(zenith * (pi / 180.0));
    size_t layers = n - 1;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < layers; ++j)
            dsdh[i * layers + j] = 0.0;
        double impact = (re + (z[layers - i] - z[0])) * sinrad;
        int id = (int)i;
        if (zenith > 90.0) {
            id = -1;
            if (impact >= re)
                for (size_t j = 1; j <= layers; ++j)
                    if (impact < (z[n - j] - z[0] + re) && impact >= (z[layers - j] - z[0] + re))
                        id = (int)j;
        }
        for (int j = 1; j <= id; ++j) {
            double sm = j == id && id == (int)i && zenith > 90.0 ? -1.0 : 1.0;
            double rj = re + (z[n - (size_t)j] - z[0]);
            double rjp1 = re + (z[layers - (size_t)j] - z[0]);
            double dhj = (z[n - (size_t)j] - z[0]) - (z[layers - (size_t)j] - z[0]);
            double ga = fmax(0.0, rj * rj - impact * impact);
            double gb = fmax(0.0, rjp1 * rjp1 - impact * impact);
            double dsj = id > (int)i && j == id ? sqrt(ga) : sqrt(ga) - sm * sqrt(gb);
            dsdh[i * layers + (size_t)j - 1] = dsj / dhj;
        }
        nid[i] = id;
    }
    return TT_OK;
}

tt_status tt_air_mass(size_t n, const double *air, const int *nid, const double *dsdh,
                      double *vertical, double *slant) {
    if (n < 2 || !air || !nid || !dsdh || !vertical || !slant)
        return TT_INVALID;
    size_t layers = n - 1;
    double accum = air[layers];
    /* TUV-x leaves the top vertical entry undefined; define it here. */
    vertical[layers] = accum;
    for (size_t i = layers; i > 0; --i) {
        accum = accum + air[i - 1];
        vertical[i - 1] = accum;
    }
    slant[layers] = dsdh[layers] * air[layers];
    for (size_t i = 1; i <= layers; ++i) {
        if (nid[i] > (int)layers)
            return TT_INVALID;
        accum = slant[layers];
        if (nid[i] < 0)
            accum = 1.e36;
        else {
            size_t single = (size_t)nid[i] < i ? (size_t)nid[i] : i;
            for (size_t j = 1; j <= single; ++j)
                accum = accum + air[layers - j] * dsdh[i * layers + j - 1];
            for (size_t j = single + 1; j <= (size_t)nid[i]; ++j)
                accum = accum + 2.0 * air[layers - j] * dsdh[i * layers + j - 1];
        }
        slant[layers - i] = accum;
    }
    return TT_OK;
}

static int leap(int year) {
    return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

int tt_day_of_year(int year, int month, int day) {
    static const int days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (year < 1950 || year > 2050 || month < 1 || month > 12 || day < 1 ||
        day > days[month - 1] + (month == 2 && leap(year)))
        return 0;
    int sum = day;
    for (int i = 1; i < month; ++i)
        sum += days[i - 1] + (i == 2 && leap(year));
    return sum;
}

static double positive_mod(double x, double y) {
    double v = fmod(x, y);
    return v < 0.0 ? v + y : v;
}

static double solar_time(int year, int day, double hour) {
    if (year < 1950 || year > 2050 || day < 1 || day > 365 + leap(year) || !isfinite(hour))
        return NAN;
    int delta = year - 1949;
    double jd = 32916.5 + (double)(delta * 365 + delta / 4 + day) + hour / 24.0;
    if (year % 100 == 0 && year % 400 != 0)
        jd = jd - 1.0;
    return jd - 51545.0;
}

double tt_earth_sun_distance(int year, int day, double hour) {
    double time = solar_time(year, day, hour);
    double d2r = (2.0 * asin(1.0)) / 180.0;
    double anomaly = positive_mod(357.528 + 0.9856003 * time, 360.0) * d2r;
    return 1.00014 - 0.01671 * cos(anomaly) - 0.00014 * cos(2.0 * anomaly);
}

double tt_solar_elevation(int year, int day, double hour, double lat, double lon) {
    if (!isfinite(lat) || !isfinite(lon) || fabs(lat) > 90.0 || fabs(lon) > 180.0)
        return NAN;
    double pi = 2.0 * asin(1.0), d2r = pi / 180.0;
    double time = solar_time(year, day, hour);
    double anomaly = positive_mod(357.528 + 0.9856003 * time, 360.0) * d2r;
    double mean_long = positive_mod(280.460 + 0.9856474 * time, 360.0);
    double obliquity = (23.439 - 0.0000004 * time) * d2r;
    double ecliptic = mean_long + 1.915 * sin(anomaly) + 0.020 * sin(2.0 * anomaly);
    ecliptic = positive_mod(ecliptic, 360.0) * d2r;
    double num = cos(obliquity) * sin(ecliptic), den = cos(ecliptic);
    double ra = atan(num / den);
    if (den < 0.0)
        ra = ra + pi;
    else if (num < 0.0)
        ra = ra + 2.0 * pi;
    double dec = asin(sin(obliquity) * sin(ecliptic));
    double gmst = positive_mod(6.697375 + 0.0657098242 * time + hour, 24.0);
    double lmst = positive_mod(gmst + lon / 15.0, 24.0);
    lmst = lmst * 15.0 * d2r;
    double ha = lmst - ra;
    if (ha < -pi)
        ha = ha + 2.0 * pi;
    else if (ha > pi)
        ha = ha - 2.0 * pi;
    return asin(sin(dec) * sin(lat * d2r) + cos(dec) * cos(lat * d2r) * cos(ha)) / d2r;
}
