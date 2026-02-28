#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>

using namespace std;

//Constants
const int n = 9;

// (L = Long double type)
static const long double p[n] = {                                               //Roots of chebisev polinomials n = 9, g = 7
    0.99999999999980993L, 676.5203681218851L, -1259.1392167224028L,
    771.32342877765313L, -176.61502916214059L, 12.507343278686905L,
    -0.13857109526572012L, 9.9843695780195716e-6L, 1.5056327351493116e-7L
};

const complex<long double> Im(0.0L, 1.0L);

const long double pi = 3.141592653589793238462643383279502884L;
const long double g = 7.0L;
const int kummerval = 1000;
const long double sigma_Bromwich = 0.0L;
const long double nu_Talbot = 1.0L;


complex<long double> Gamma(complex<long double> z) {                    //Gamma function: https://en.wikipedia.org/wiki/Lanczos_approximation
    if (real(z) < 0.5L) {
        return pi / (sin(pi * z) * Gamma(1.0L - z));                //Reflection formula
    }
    else {
        z -= 1.0L;
        complex<long double> x = p[0];                  //Stirling's fourmula when z is big Lanczos's approximation
        for (int i = 1; i < n; i++) {
            x += p[i] / (z + (long double)i);
        }
        complex<long double> t = z + g + 0.5L;
        return sqrt(2.0L * pi) * pow(t, (z + 0.5L)) * exp(-t) * x;
    }
}
// https://en.wikipedia.org/wiki/Confluent_hypergeometric_function#Asymptotic_behavior: Taylor expansion
complex<long double> Kummer(complex<long double> a, complex<long double> b, long double z) {
    // Kummer transformation for stabiliy
    // if z is negative we transform it to positive
    if (z < 0.0L) {
        return exp(z) * Kummer(b - a, b, -z);
    }

    // Taylor expansion, stable because z is bigger than 0
    complex<long double> sum(1.0L, 0.0L);
    complex<long double> term(1.0L, 0.0L);

    for(int i = 1; i < kummerval; i++) {
        term = term * ((a + (long double)i - 1.0L) / (b + (long double)i - 1.0L)) * (z / (long double)i);
        sum = sum + term;

            // if term is smaall we can end it

    }
    return sum;
}

complex<long double> Laplace(complex<long double> lambda, long double T, long double K, long double nu) {
    complex<long double> gamma_val = sqrt(2.0L * lambda + nu * nu);
    complex<long double> alpha_val = (gamma_val + nu) / 2.0L;
    complex<long double> beta_val = (gamma_val - nu) / 2.0L;



    return (pow(2.0L * K, (1.0L - beta_val))) / (2.0L * lambda * (beta_val - 1.0L)) * (Gamma(alpha_val + 1.0L) / Gamma(alpha_val + beta_val + 1.0L)) * Kummer(beta_val - 1.0L, alpha_val + beta_val + 1.0L, -1.0L / (2.0L * K));
}



complex<long double> TalbotCurve(long double mu_Bromwich, long double theta) {
    return sigma_Bromwich + mu_Bromwich * (theta / tan(theta) + nu_Talbot * Im * theta); //by definition, talbot curve
}

complex<long double> derivTalbotCurve(long double mu_Bromwich, long double theta) {
    return mu_Bromwich * (1.0L / tan(theta) - theta / (sin(theta) * sin(theta)) + nu_Talbot * Im); //derivative then L'Hopital's rule 0/0 lim_theta \rightarrow 0
}

complex<long double> IntglTalbot(int M, long double T, long double K, long double r, long double sigma) {

    long double nu = (2.0L * r / (sigma * sigma)) - 1.0L;

    long double mu_Bromwich = 2.0L * M / (5.0L * T);
    long double incr_intgr = 2.0L * pi / M; // from -pi to pi the integral is, so the rectangles will have bases this long.

    long double theta_k = 0.0L;
    complex<long double> sumint(0.0L, 0.0L);

    for (int i = 0; i < M; i++) {
        // Midpoint rule, so the integral is calculated using the midpoint on the x axis using that value to calculate rectangles height
        theta_k = -pi + incr_intgr * ((long double)i + 0.5L); //- 0.5 so it will be on the middle points

        complex<long double> talbotkont = TalbotCurve(mu_Bromwich, theta_k);
        complex<long double> d_talbotkont = derivTalbotCurve(mu_Bromwich, theta_k);

        sumint += exp(talbotkont * T) * Laplace(talbotkont, T, K, nu) * d_talbotkont;   //integral expressed as sum of rectangles
    }

    return (sumint * incr_intgr) / (2.0L * pi * Im);
}



long double OptionPrice(long double S_0, long double T, long double K, long double sigma, long double r) {
    int M = 200;


    complex<long double> result = IntglTalbot(M,  T * (sigma * sigma) / 4.0L, K * T * (sigma * sigma) / (4.0L * S_0), r, sigma);

    return exp(-r * T) * (S_0 / T) * (4.0L / (sigma * sigma)) * real(result);
}

// --- 5. MAIN ---

int main() {
//one option's price
    long double S_0 = 100.0L;
    long double K = 100.0L;
    long double T = 1.0L;
    long double sigma = 0.11L;
    long double r = 0.05L;






    long double price = OptionPrice(S_0, T, K, sigma, r);

    cout << price << endl;

//Ai generated code : Random number generator, using my instructions for ai training, dangerous when sigma^2*T < 0.01, and export to csv file

/*
    const int NUM_SAMPLES = 10000;   // Hány adatot szeretnél?
    ofstream file("training_data.csv"); // Kimeneti fájl

    // CSV fejléc (S0, K, T, sigma, r, price)
    file << "S0,K,T,sigma,r,price" << endl;

    // 2. MODERN VÉLETLENSZÁM GENERÁTOR INICIALIZÁLÁSA
    random_device rd;  // "Hardveres" mag a véletlenszámokhoz
    mt19937 gen(rd()); // Mersenne Twister generátor (nagyon jó minõségû)

    // 3. INTERVALLUMOK BEÁLLÍTÁSA (A DOKUMENTUM ALAPJÁN)
    // long double típust használunk a pontosság miatt
    uniform_real_distribution<long double> dist_S0(80.0, 120.0);    // Árfolyam: 80 - 120
    uniform_real_distribution<long double> dist_K(80.0, 120.0);     // Kötési ár: 80 - 120

    // Figyelem: A T-t ne engedjük 0.25 (3 hónap) alá a stabilitás miatt
    uniform_real_distribution<long double> dist_T(0.25, 2.0);       // Lejárat: 0.25 év - 2 év

    // Figyelem: A volatilitást ne engedjük 0.1 alá a "veszélyzóna" miatt
    uniform_real_distribution<long double> dist_sigma(0.12, 0.50);  // Volatilitás: 12% - 50%

    uniform_real_distribution<long double> dist_r(0.01, 0.10);      // Kamat: 1% - 10%

    cout << "Adatgeneralas inditasa..." << endl;
    int valid_samples = 0;

    // 4. GENERÁLÓ CIKLUS (While loop, hogy biztosan meglegyen a 10.000 jó adat)
    while (valid_samples < NUM_SAMPLES) {

        // Paraméterek sorsolása
        long double S0 = dist_S0(gen);
        long double K = dist_K(gen);
        long double T = dist_T(gen);
        long double sigma = dist_sigma(gen);
        long double r = dist_r(gen);

        // --- KRITIKUS SZÛRÕ A DOKUMENTUM ALAPJÁN ---
        // A tanulmány szerint, ha sigma^2 * T < 0.01, az inverzió instabil lehet. [cite: 400]
        if (sigma * sigma * T < 0.01) {
            continue; // Ezt a kört eldobjuk, generálunk újat
        }

        // Árazás
        long double price = OptionPrice(S0, T, K, sigma, r);

        // Ellenõrzés: Ha az eredmény "Not a Number" vagy végtelen, eldobjuk
        if (isnan(price) || isinf(price) || price < 0.0L || price > S0) {
            continue;
        }

        // Ha minden oké, kiírjuk a fájlba
        file << S0 << "," << K << "," << T << "," << sigma << "," << r << "," << price << endl;

        valid_samples++;

        // Progress jelzõ (hogy lásd, hol tart)
        if (valid_samples % 1000 == 0) {
            cout << valid_samples << " / " << NUM_SAMPLES << " kesz..." << endl;
        }
    }

    file.close();
    cout << "Kesz! Az adatok a training_data.csv fajlban." << endl;
*/

    return 0;
}
