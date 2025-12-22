#include <iostream>
#include <iomanip>

double utc_to_mjd_tdb(int year, int month, int day, int hour, int min, double sec) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12*a - 3;
    double jd = day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
    jd += (hour - 12)/24.0 + min/1440.0 + sec/86400.0;
    double mjd = jd - 2400000.5;
    return mjd + 69.184 / 86400.0;  // TT = UTC + 69.184s approx
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Epoca elementi: MJD 61000.000000\n";
    for (int hour = 0; hour <= 24; hour += 4) {
        double mjd_tdb = utc_to_mjd_tdb(2025, 11, 28, hour, 0, 0.0);
        double delta = mjd_tdb - 61000.0;
        std::cout << "Ora " << std::setw(2) << hour << ":00 → MJD " << mjd_tdb 
                  << " (Δ = " << delta << " giorni)\n";
    }
    return 0;
}
