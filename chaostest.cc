#include <stdlib.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
// adapted from the algorithm of https://hal.archives-ouvertes.fr/hal-02388470/document
// and remarks from
// https://www.maths.usyd.edu.au/u/gottwald/preprints/testforchaos_MPI.pdf
// the notation are the same as in those articles except phi that became data
// https://math.stackexchange.com/questions/2603548/least-absolute-deviation-lad-line-fitting-regression
// define Pi and the number of different value of c
// https://www.cec.uchile.cl/cinetica/pcordero/MC_libros/NumericalRecipesinC.pdf
#define M_PI 3.14159265358979323846 /* pi */
#define NC 100
#define THRESOLD 0.00000001  // used for linear regression
// read a file where there is no empty line and each line has this form :time data\n

void read_file(std::string const& filename, std::vector<double>& time_file,
               std::vector<double>& data_file) {
    std::string line;
    std::ifstream file(filename);
    std::string::size_type pos;
    int i = 0;
    if (file.is_open()) {
        while (std::getline(file, line)) {
            if (i > 1000) {
                time_file.push_back((double)(i + 1));
                std::stod(line, &pos);
                try {
                    data_file.push_back((std::stod(line.substr(pos))));
                } catch (const std::out_of_range& e) {
                    data_file.push_back(0);
                }
            }
            i++;
        }
        file.close();
    } else {
        std::cout << "impossible to open file" << std::endl;
        exit(2);
    }
}

double average(std::vector<double> const& v) {
    if (v.empty()) {
        return 0;
    }

    // double const count = static_cast<double>(v.size());
    return std::accumulate(v.begin(), v.end(), 0) / (double)v.size();
}

double covariance(
    std::vector<double> const& a, std::vector<double> const& b,
    const std::size_t n) {  // we are not computing all value of M so if a is M and b is
                            // the time the vector will not have the same size.
    double meana = average(a);
    double meanb = average(b);
    double res = 0;
    for (std::size_t i = 1; i < n; i++) {
        res += (a[i] - meana) * (b[i] - meanb);
    }
    return res / n;
}

double compute_Kn_correlation(std::vector<double> const& a, std::vector<double> const& b,
                              const std::size_t n0) {
    return std::abs(covariance(a, b, n0) /
                    std::sqrt(covariance(a, a, n0) * covariance(b, b, n0)));
}

// the regression method use the least deviation as said in the first article. The choosen
// algorithme is the iteratively reweighted least square to understand this function you
// can do all calculation (since p=1 for IRLS) or read the joined pdf if there is one.
// regression :ax+b time need to start at a value >0
void compute_Kn_regression(std::vector<double> const& log_time,
                           std::vector<double> const& log_Dtilde, const std::size_t n0,
                           double& ar) {
    double a = 1.0, b = 0.0;  // init
    double a_old = 0.0,
           b_old =
               0.0;  // keep value of the precedent iteration (too check the convergence)
    double w_i = 1.0;             // weight
    double threshold_0 = 0.0001;  // to avoid dividing by zeros
    double sum_wi, sum_wi_phii, sum_wi_phii_ti, sum_wi_ti, sum_wi_ti2;
    while (std::pow(a - a_old, 2) + std::pow(b - b_old, 2) > THRESOLD) {
        a_old = a;
        b_old = b;
        sum_wi = 0;
        sum_wi_phii = 0;
        sum_wi_phii_ti = 0;
        sum_wi_ti = 0;
        sum_wi_ti2 = 0;
        for (std::size_t i = 0; i < n0; i++) {
            w_i = 1.0 /
                  (std::max(std::abs(log_Dtilde[i] - a * log_time[i] - b), threshold_0));
            sum_wi += w_i;
            sum_wi_phii += w_i * log_Dtilde[i];
            sum_wi_phii_ti += w_i * log_Dtilde[i] * log_time[i];
            sum_wi_ti += w_i * log_time[i];
            sum_wi_ti2 += w_i * std::pow(log_time[i], 2);

        }
        a = (1 / (sum_wi_ti2 * sum_wi - sum_wi_ti * sum_wi_ti)) *
            (sum_wi * sum_wi_phii_ti - sum_wi_ti * sum_wi_phii);

        b = (1 / (sum_wi_ti2 * sum_wi - sum_wi_ti * sum_wi_ti)) *
            (sum_wi_ti2 * sum_wi_phii - sum_wi_ti * sum_wi_phii_ti);
    }
    ar = a;
}


void compute_pq(double& c, std::vector<double>& data, std::vector<double>& p,
                std::vector<double>& q, std::size_t n) {
    p[0] = data[0] * std::cos(c);
    q[0] = data[0] * std::sin(c);
    for (std::size_t i = 1; i < n; i++) {
        p[i] = p[i - 1] + data[i] * std::cos(c * (i + 1));
        q[i] = q[i - 1] + data[i] * std::cos(c * (i + 1));
    }
}

void compute_M(std::vector<double>& M, const std::vector<double>& p,
               const std::vector<double>& q, const std::size_t n, const std::size_t n0) {
    for (std::size_t i = 0; i < n0; i++) {
        M[i] = 0.0;
        for (std::size_t j = 0; j + i < n; j++) {
            M[i] += std::pow(p[j + i] - p[j], 2) + std::pow(q[j + i] - q[j], 2);
        }
        M[i] /= (n - i);
    }
}

void compute_D(double& Edata2, std::vector<double>& D, const std::vector<double>& p,
               const std::vector<double>& q, const std::size_t& n, const std::size_t n0,
               const double& c) {
    for (std::size_t i = 0; i < n0; i++) {
        D[i] = 0.0;
        for (std::size_t j = 0; j + i < n; j++) {
            D[i] += std::pow(p[j + i] - p[j], 2) + std::pow(q[j + i] - q[j], 2);
        }
        D[i] /= (n - i);
        D[i] -= Edata2 * (1 - std::cos((i + 1) * c)) / (1 - std::cos(c));
    }
}

// replace D by Dtile (we do not need D anymore so we use the space available to store
// Dtilde)
void compute_Dtilde(std::vector<double>& D, const double& a) {
    double min = std::abs(*std::min_element(D.begin(), D.end()));
    std::size_t n = D.size();
    if (min == 0) min += THRESOLD;  // used for linear regression
    ;
    for (std::size_t i = 0; i < n; i++) {
        D[i] += a * min;
    }
}

double median(std::vector<double> v) {
    std::size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

double correlation_method(std::vector<double>& time, std::vector<double>& data,
                          std::vector<double>& pc, std::vector<double>& qc,
                          std::vector<double>& D, std::vector<double>& K, std::size_t& n,
                          std::size_t& N0) {
    std::mt19937_64 rng;  // Mersenne Twister generator
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);
    // initialize a uniform distribution between
    std::uniform_real_distribution<double> unif(0.1, 2 * M_PI - 0.1);
    double delta = M_PI / (NC + 1);
    double c = delta;
    double Edata2 = std::pow(average(data), 2);
    /////////Loop computation/////////
    for (int ccnt = 0; ccnt < NC; ccnt++) {  // the loop index is ccnt for c counter
        // std::cout << ccnt+1 << "/" << NC << std::endl;
        // c += delta;
        c = unif(rng);
        compute_pq(c, data, pc, qc, n);
        // compute_M(D, pc, qc, n, N0);
        compute_D(Edata2, D, pc, qc, n, N0, c);
        // double Kn=compute_Kn_correlation(time,D,N0);
        K[ccnt] = compute_Kn_correlation(time, D, N0);
        // std::cout << "coorr " << K[ccnt] << std::endl;
    }
    return median(K);
}

double regression_method(std::vector<double>& time, std::vector<double>& data,
                         std::vector<double>& pc, std::vector<double>& qc,
                         std::vector<double>& D, std::vector<double>& K, std::size_t& n,
                         std::size_t& N0) {
    std::mt19937_64 rng;  // Mersenne Twister generator
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);
    // initialize a uniform distribution between
    std::uniform_real_distribution<double> unif(0.1, 2 * M_PI - 0.1);
    // double delta = (3 * M_PI / 5) / (NC + 1);
    // double c = M_PI / 5;
    double delta = M_PI / (NC + 1);
    double c = delta;
    double Edata2 = std::pow(average(data), 2);
    std::vector<double> log_time(n);
    std::vector<double> log_dtilde(N0);
    for (size_t i = 0; i < n; i++) {
        log_time[i] = std::log(time[i]);
    }
    /////////Loop computation/////////
    for (int ccnt = 0; ccnt < NC; ccnt++) {  // the loop index is ccnt for c counter
        // c += delta;
        c = unif(rng);
        // std::cout<< "c: " << c<< std::endl;
        compute_pq(c, data, pc, qc, n);
        compute_D(Edata2, D, pc, qc, n, N0, c);
        compute_Dtilde(D, 1.1);  // D is becoming dtilde
        for (size_t i = 0; i < N0; i++) {
            log_dtilde[i] = std::log(D[i]);
        }
        compute_Kn_regression(log_time, log_dtilde, N0, K[ccnt]);
        if (std::isnan(K[ccnt])) {
            std::cout << "NAN was found" << NC << std::endl;
        }
        // std::cout << ccnt+1 << "/" << NC << std::endl;
    }
    return median(K);
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "chaostest file" << std::endl;
        exit(1);
    }
    std::string path(argv[1]);
    /////////defining usefull variable/////////
    std::vector<double> time;
    std::vector<double> data;
    read_file(path, time, data);
    std::size_t n = data.size();
    std::vector<double> pc(n);
    std::vector<double> qc(n);
    std::vector<double> K(NC);
    std::size_t N0 = n / 10;
    std::vector<double> D(N0);
    double res = correlation_method(time, data, pc, qc, D, K, n, N0);
    std::cout << "score by correlation:" << res << std::endl;
    res = regression_method(time, data, pc, qc, D, K, n, N0);
    std::cout << "score by regression:" << res << std::endl;
    return 0;
}

