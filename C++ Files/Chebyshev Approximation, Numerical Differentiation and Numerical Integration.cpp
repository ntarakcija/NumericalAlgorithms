#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iomanip>

double const pi(atan(1) * 4);

class ChebyshevApproximation {
    int n, m;
    double xmin, xmax;
    std::vector<double> c;
    public:
    template <typename FunTip>
    ChebyshevApproximation(FunTip f, double xmin, double xmax, int n);
    void set_m(int m);
    void trunc(double eps);
    double operator()(double x) const;
    double derivative(double x) const;
    ChebyshevApproximation derivative() const;
    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

template <typename FunTip>
ChebyshevApproximation::ChebyshevApproximation(FunTip f, double xmin, double xmax, int n) : n(n), m(n), xmin(xmin), xmax(xmax) {
    if(n < 1 || xmin >= xmax) throw std::domain_error("Bad parameters");
    
    c.resize(n + 1);
    std::vector<double> w(n + 2);
    std::vector<double> v(n + 1);
    for(int i(0); i <= n + 1; i++)
        w[i] = cos(pi * i / (2 * n + 2));
    for(int i(0); i <= n / 2; i++)
        v[i] = f((xmin + xmax + (xmax - xmin) * w[2 * i + 1]) / 2);
    for(int i(n / 2 + 1); i <= n; i++)
        v[i] = f((xmin + xmax + (xmax - xmin) * w[2 * n + 1 - 2 * i]) / 2);
    for(int k(0); k <= n; k++) {
        double s(0);
        for(int i(0); i <= n; i++) {
            int p((k * (2 * i + 1)) % (4 * n + 4));
            if(p >= 2 * n + 2)
                p = 4 * n + 4 - p;
            if(p >= n + 1)
                s -= v[i] * w[2 * n + 2 - p];
            else
                s += v[i] * w[p];
        }
        c[k] = 2 * s / (n + 1);
    }
}

inline void ChebyshevApproximation::set_m(int m) {
    if(m <= 1 || m > n) throw std::domain_error("Bad order");
    else this->m = m;
}

inline void ChebyshevApproximation::trunc(double eps) {
    if(eps < 0) throw std::domain_error("Bad tolerance");
    
    for(int i(c.size() - 1); i >= 0; i--)
        if(std::abs(c[i]) < eps) m--;
        else break;
        
    if(m < 1) throw std::domain_error("Bad tolerance");
}

double ChebyshevApproximation::operator()(double x) const {
    if(x < xmin || x > xmax) throw std::domain_error("Bad argument");
    
    double t((2 * x - xmin - xmax) / (xmax - xmin));
    double p(1);
    double q(t);
    double s(c[0] / 2 + c[1] * t);
    for(int k(2); k <= m; k++) {
        double r(2 * t * q - p);
        s += c[k] * r;
        p = q;
        q = r;
    }
    return s;
}

double ChebyshevApproximation::derivative(double x) const {
    if(x < xmin || x > xmax) throw std::domain_error("Bad argument");
    
    double t((2 * x - xmin - xmax) / (xmax - xmin));
    double p(1);
    double q(4 * t);
    double s(c[1] + 4 * c[2] * t);
    for(int k(3); k <= m; k++) {
        double r(k * (2 * t * q / (k - 1) - p / (k - 2)));
        s += c[k] * r;
        p = q;
        q = r;
    }
    return 2 * s / (xmax - xmin);
}

ChebyshevApproximation ChebyshevApproximation::derivative() const {
    std::vector<double> cprim(m + 1);
    double mi(4 / (xmax - xmin));
    cprim[m - 1] = mi * m * c[m];
    cprim[m - 2] = mi * (m - 1) * c[m - 1];
    for(int k(m - 3); k >= 0; k--)
        cprim[k] = cprim[k + 2] + mi * (k + 1) * c[k + 1];
    ChebyshevApproximation ca(*this);
    ca.c = cprim;
    return ca;
}

ChebyshevApproximation ChebyshevApproximation::antiderivative() const {
    std::vector<double> czvijezda(m + 3);
    czvijezda[m + 2] = 0;
    czvijezda[m + 1] = 0;
    for(int k(m); k >= 0; k--)
        czvijezda[k] = ((xmax - xmin) / (4 * k)) * (c[k - 1] - c[k + 1]);
    ChebyshevApproximation ca(*this);
    ca.c = czvijezda;
    return ca;
}

double ChebyshevApproximation::integrate(double a, double b) const {
    if(a < xmin || b < xmin || a > xmax || b > xmax) throw std::domain_error("Bad interval");

    double mi((b - a) / 2);
    double suma(0);
    for(int k(1); k <= (m + 1) / 2; k++)
        suma += (2 * c[2 * k]) / (1 - 4 * k * k);
    if(a < b) return mi * c[0] + mi * suma;
    else return -(mi * c[0] + mi * suma);
}

double ChebyshevApproximation::integrate() const {
    double mi((xmax - xmin) / 2);
    double suma(0);
    for(int k(1); k <= (m + 1) / 2; k++)
        suma += (2 * c[2 * k]) / (1 - 4 * k * k);
    if(xmin < xmax) return mi * c[0] + mi * suma;
    else return -(mi * c[0] + mi * suma);
}


template <typename FunTip>
std::pair<double, bool> RombergIntegration(FunTip f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 50) {
    if(eps < 0 || nmin < 0 || nmax < 0 || nmin > nmax) throw std::domain_error("Bad parameter");
    
    std::pair<double, bool> rjesenje;
    int Mmax(50);
    int N(2);
    double h((b - a) / N);
    double s((f(a) + f(b)) / 2);
    double Iold(s);
    std::vector<double> I(Mmax + 1);
    for(int i(1); i <= Mmax; i++) {
        for(int j(1); j <= N / 2; j++)
            s += f(a + (2 * j - 1) * h);
        I[i] = h * s;
        double p(4);
        for(int k(i - 1); k >= 1; k--) {
            I[k] = (p * I[k + 1] - I[k]) / (p - 1);
            p *= 4;
        }
        if(std::abs(I[1] - Iold) <= eps && N > nmin) {
            if(a > b) rjesenje.first = -I[1];
            else rjesenje.first = I[1];
            rjesenje.second = true;
            return rjesenje;
        }
        Iold = I[1];
        h /= 2;
        N *= 2;
        if(N >= nmax) break;
    }
    
    if(a > b) rjesenje.first = -I[1];
    else rjesenje.first = I[1];
    rjesenje.second = false;
    return rjesenje;
}

template <typename FunTip>
std::pair<double, bool> TanhSinhIntegration(FunTip f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 20, double range = 3.5) {
    if(eps < 0 || nmin < 0 || nmax < 0 || nmin > nmax || range < 0) throw std::domain_error("Bad parameter");
    
    // Modifikacije funkcije f
    auto fun = [f, a, b](double x) {
        double temp = (cosh(x) / (cosh((pi / 2) * sinh(x)) * cosh((pi / 2) * sinh(x)))) * f((b + a) / 2 + ((b - a) / 2) * tanh((pi / 2) * sinh(x)));
        if(std::isfinite(temp)) return temp;
        else return 0.;
    };
    
    std::pair<double, bool> rjesenje;
    double a1(-range), b1(range);
    int N(2);
    double h((b1 - a1) / N);
    double s((fun(a1) + fun(b1)) / 2);
    double Iold(s);
    double I;
    while(N < nmax) {
        for(int i(1); i <= N / 2; i++)
            s += fun(a1 + (2 * i - 1) * h);
        I = h * s;
        if(N > nmin && std::abs(I - Iold) <= eps) {
            rjesenje.first = I;
            rjesenje.second = true;
            break;
        }
        Iold = I;
        N *= 2;
        h /= 2;
    }
    rjesenje.first = I;
    if(N > nmax) rjesenje.second = false;
    
    rjesenje.first *= ((b - a) * pi / 4);
    if(a > b) rjesenje.first *= -1; 
    return rjesenje;
}

bool presaoDubinu(false);
template <typename FunTip>
double AdaptiveAux(FunTip f, double a, double b, double eps, double f1, double f2, double f3, int R) {
    double c((a + b) / 2);
    double I1((b - a) * (f1 + 4 * f3 + f2) / 6);
    double f4(f((a + c) / 2));
    double f5(f((c + b) / 2));
    double I2((b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12);
    if(R <= 0 || std::abs(I1 - I2) <= eps) {
        if(R <= 0) presaoDubinu = true;
        return I2;
    }
    return AdaptiveAux(f, a, c, eps, f1, f3, f4, R - 1) + AdaptiveAux(f, c, b, eps, f3, f2, f5, R - 1);
}

template <typename FunTip>
std::pair<double, bool> AdaptiveIntegration(FunTip f, double a, double b, double eps = 1e-10, int maxdepth = 30, int nmin = 1) {
    if(eps < 0 || maxdepth < 0 || nmin < 0) throw std::domain_error("Bad parameter");
    
    // Modifikacija funkcije f za slučaj kada je vrijednost beskonačno
    auto fun = [f](double x) {
        double temp = f(x);
        if(std::isfinite(temp)) return temp;
        else return 0.;
    };
    
    double s(0);
    double h((b - a) / nmin);
    std::pair<double, bool> rjesenje;
    for(int i(1); i <= nmin; i++) {
        presaoDubinu = false;
        s += AdaptiveAux(fun, a, a + h, eps, fun(a), fun(a + h), fun(a + h / 2), maxdepth); 
        if(presaoDubinu) rjesenje.second = false;
        else rjesenje.second = true;
        a += h;
    }
    
    rjesenje.first = s;
    if(a > b) rjesenje.first *= -1;
    return rjesenje;
}


int main () {
    // Testiranje klase Chebyshev Approximation
    auto sinus = [](double x) { return std::sin(x); };
    ChebyshevApproximation ca1(sinus, 0, pi, 1000);
    std::cout << std::setw(10) << std::left << "Tacka";
    std::cout << std::setw(15) << std::left << "Sinus";
    std::cout << std::setw(15) << std::left << "Aproksimacija";
    std::cout << std::setw(15) << std::left << "Broj uzoraka: 1000" << std::endl;
    std::cout << std::setw(10) << std::left << "0";
    std::cout << std::setw(15) << std::left << sinus(0);
    std::cout << std::setw(15) << std::left << ca1(0) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/3";
    std::cout << std::setw(15) << std::left << sinus(pi / 3);
    std::cout << std::setw(15) << std::left << ca1(pi / 3) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/4";
    std::cout << std::setw(15) << std::left << sinus(pi / 4);
    std::cout << std::setw(15) << std::left << ca1(pi / 4) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/6";
    std::cout << std::setw(15) << std::left << sinus(pi / 6);
    std::cout << std::setw(15) << std::left << ca1(pi / 6) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/2";
    std::cout << std::setw(15) << std::left << sinus(pi / 2);
    std::cout << std::setw(15) << std::left << ca1(pi / 2) << std::endl;
    std::cout << std::setw(10) << std::left << "pi";
    std::cout << std::setw(15) << std::left << sinus(pi);
    std::cout << std::setw(15) << std::left << ca1(pi) << std::endl;
    std::cout << std::endl;
    
    ca1.trunc(1e-5); // Odbacivanje nekih koeficijenata - trebala bi se smanjiti tacnost
    std::cout << std::setw(10) << std::left << "Tacka";
    std::cout << std::setw(15) << std::left << "Sinus";
    std::cout << std::setw(15) << std::left << "Aproksimacija";
    std::cout << std::setw(15) << std::left << "Broj uzoraka: trunc" << std::endl;
    std::cout << std::setw(10) << std::left << "0";
    std::cout << std::setw(15) << std::left << sinus(0);
    std::cout << std::setw(15) << std::left << ca1(0) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/3";
    std::cout << std::setw(15) << std::left << sinus(pi / 3);
    std::cout << std::setw(15) << std::left << ca1(pi / 3) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/4";
    std::cout << std::setw(15) << std::left << sinus(pi / 4);
    std::cout << std::setw(15) << std::left << ca1(pi / 4) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/6";
    std::cout << std::setw(15) << std::left << sinus(pi / 6);
    std::cout << std::setw(15) << std::left << ca1(pi / 6) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/2";
    std::cout << std::setw(15) << std::left << sinus(pi / 2);
    std::cout << std::setw(15) << std::left << ca1(pi / 2) << std::endl;
    std::cout << std::setw(10) << std::left << "pi";
    std::cout << std::setw(15) << std::left << sinus(pi);
    std::cout << std::setw(15) << std::left << ca1(pi) << std::endl;
    std::cout << std::endl;
    
    ca1.set_m(5); // Smanjivanje uzoraka - trebala bi se smanjiti tačnost
    std::cout << std::setw(10) << std::left << "Tacka";
    std::cout << std::setw(15) << std::left << "Sinus";
    std::cout << std::setw(15) << std::left << "Aproksimacija";
    std::cout << std::setw(15) << std::left << "Broj uzoraka: 5" << std::endl;
    std::cout << std::setw(10) << std::left << "0";
    std::cout << std::setw(15) << std::left << sinus(0);
    std::cout << std::setw(15) << std::left << ca1(0) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/3";
    std::cout << std::setw(15) << std::left << sinus(pi / 3);
    std::cout << std::setw(15) << std::left << ca1(pi / 3) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/4";
    std::cout << std::setw(15) << std::left << sinus(pi / 4);
    std::cout << std::setw(15) << std::left << ca1(pi / 4) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/6";
    std::cout << std::setw(15) << std::left << sinus(pi / 6);
    std::cout << std::setw(15) << std::left << ca1(pi / 6) << std::endl;
    std::cout << std::setw(10) << std::left << "pi/2";
    std::cout << std::setw(15) << std::left << sinus(pi / 2);
    std::cout << std::setw(15) << std::left << ca1(pi / 2) << std::endl;
    std::cout << std::setw(10) << std::left << "pi";
    std::cout << std::setw(15) << std::left << sinus(pi);
    std::cout << std::setw(15) << std::left << ca1(pi) << std::endl;
    std::cout << std::endl;
    
    auto kubna = [](double x) { return x * x * x; };
    ChebyshevApproximation ca2(kubna, -100, 100, 1000);
    std::cout << std::setw(10) << std::left << "Tacka";
    std::cout << std::setw(15) << std::left << "Kubna";
    std::cout << std::setw(15) << std::left << "Aproksimacija";
    std::cout << std::setw(15) << std::left << "Broj uzoraka: 1000" << std::endl;
    std::cout << std::setw(10) << std::left << "-17";
    std::cout << std::setw(15) << std::left << kubna(-17);
    std::cout << std::setw(15) << std::left << ca2(-17) << std::endl;
    std::cout << std::setw(10) << std::left << "-5";
    std::cout << std::setw(15) << std::left << kubna(-5);
    std::cout << std::setw(15) << std::left << ca2(-5) << std::endl;
    std::cout << std::setw(10) << std::left << "0";
    std::cout << std::setw(15) << std::left << kubna(0);
    std::cout << std::setw(15) << std::left << ca2(0) << std::endl;
    std::cout << std::setw(10) << std::left << "5";
    std::cout << std::setw(15) << std::left << kubna(5);
    std::cout << std::setw(15) << std::left << ca2(5) << std::endl;
    std::cout << std::setw(10) << std::left << "10";
    std::cout << std::setw(15) << std::left << kubna(10);
    std::cout << std::setw(15) << std::left << ca2(10) << std::endl;
    std::cout << std::setw(10) << std::left << "17";
    std::cout << std::setw(15) << std::left << kubna(17);
    std::cout << std::setw(15) << std::left << ca2(17) << std::endl;
    std::cout << std::endl;
    
    ca2.trunc(1e-1); // Odbacivanje nekih koeficijenata - trebala bi se smanjiti tacnost
    std::cout << std::setw(10) << std::left << "Tacka";
    std::cout << std::setw(15) << std::left << "Kubna";
    std::cout << std::setw(15) << std::left << "Aproksimacija";
    std::cout << std::setw(15) << std::left << "Broj uzoraka: trunc" << std::endl;
    std::cout << std::setw(10) << std::left << "-17";
    std::cout << std::setw(15) << std::left << kubna(-17);
    std::cout << std::setw(15) << std::left << ca2(-17) << std::endl;
    std::cout << std::setw(10) << std::left << "-5";
    std::cout << std::setw(15) << std::left << kubna(-5);
    std::cout << std::setw(15) << std::left << ca2(-5) << std::endl;
    std::cout << std::setw(10) << std::left << "0";
    std::cout << std::setw(15) << std::left << kubna(0);
    std::cout << std::setw(15) << std::left << ca2(0) << std::endl;
    std::cout << std::setw(10) << std::left << "5";
    std::cout << std::setw(15) << std::left << kubna(5);
    std::cout << std::setw(15) << std::left << ca2(5) << std::endl;
    std::cout << std::setw(10) << std::left << "10";
    std::cout << std::setw(15) << std::left << kubna(10);
    std::cout << std::setw(15) << std::left << ca2(10) << std::endl;
    std::cout << std::setw(10) << std::left << "17";
    std::cout << std::setw(15) << std::left << kubna(17);
    std::cout << std::setw(15) << std::left << ca2(17) << std::endl;
    std::cout << std::endl;
    
    ca2.set_m(5); // Smanjivanje uzoraka - trebala bi se smanjiti tačnost
    std::cout << std::setw(10) << std::left << "Tacka";
    std::cout << std::setw(15) << std::left << "Kubna";
    std::cout << std::setw(15) << std::left << "Aproksimacija";
    std::cout << std::setw(15) << std::left << "Broj uzoraka: 5" << std::endl;
    std::cout << std::setw(10) << std::left << "-17";
    std::cout << std::setw(15) << std::left << kubna(-17);
    std::cout << std::setw(15) << std::left << ca2(-17) << std::endl;
    std::cout << std::setw(10) << std::left << "-5";
    std::cout << std::setw(15) << std::left << kubna(-5);
    std::cout << std::setw(15) << std::left << ca2(-5) << std::endl;
    std::cout << std::setw(10) << std::left << "0";
    std::cout << std::setw(15) << std::left << kubna(0);
    std::cout << std::setw(15) << std::left << ca2(0) << std::endl;
    std::cout << std::setw(10) << std::left << "5";
    std::cout << std::setw(15) << std::left << kubna(5);
    std::cout << std::setw(15) << std::left << ca2(5) << std::endl;
    std::cout << std::setw(10) << std::left << "10";
    std::cout << std::setw(15) << std::left << kubna(10);
    std::cout << std::setw(15) << std::left << ca2(10) << std::endl;
    std::cout << std::setw(10) << std::left << "17";
    std::cout << std::setw(15) << std::left << kubna(17);
    std::cout << std::setw(15) << std::left << ca2(17) << std::endl;
    std::cout << std::endl;
    
    ca1.set_m(1000);
    std::cout << ca1.derivative(pi/2) << std::endl; // sin(x)' = cos(x) -> cos(pi/2) = 0
    std::cout << ca1.derivative(0) << std::endl; // sin(x)' = cos(x) -> cos(0) = 1
    std::cout << ca1.derivative(pi/4) << std::endl; // sin(x)' = cos(x) -> cos(pi/4) = sqrt(2)/2
    std::cout << std::endl;
    
    std::cout << ca1.derivative()(pi/2) << std::endl;
    std::cout << ca1.derivative()(0) << std::endl;
    std::cout << ca1.derivative()(pi/4) << std::endl;
    std::cout << std::endl;
    
    std::cout << ca1.derivative().derivative(pi/2) << std::endl; // sin(x)'' = cos(x)' = -sin(x) -> -sin(pi/2) = -1
    std::cout << ca1.derivative().derivative(0) << std::endl; // sin(x)'' = cos(x)' = -sin(x) -> -sin(0) = 0
    std::cout << ca1.derivative().derivative(pi/4) << std::endl; // sin(x)'' = cos(x)' = -sin(x) -> -sin(pi/4) = -sqrt(2)/2
    std::cout << std::endl;
    
    std::cout << ca1.integrate(0, pi/2) << std::endl; // 1
    std::cout << ca1.integrate(pi/2, pi) << std::endl; // 1
    std::cout << ca1.integrate(0, pi) << std::endl; // 2
    std::cout << ca1.integrate() << std::endl; // 2
    std::cout << std::endl;
    
    ca2.set_m(1000);
    std::cout << ca2.derivative(0) << std::endl; // x^3' = 2x^2 -> 3x^2(0) = 0
    std::cout << ca2.derivative(-10) << std::endl; // x^3' = 3x^2 -> 3x^2(10) = -300
    std::cout << ca2.derivative(10) << std::endl; // x^3' = 3x^2 -> 3x^2(10) = 300
    std::cout << std::endl;
    
    std::cout << ca2.derivative()(0) << std::endl;
    std::cout << ca2.derivative()(-10) << std::endl;
    std::cout << ca2.derivative()(10) << std::endl;
    std::cout << std::endl;
    
    std::cout << ca2.derivative().derivative(0) << std::endl; // x^3'' = 3x^2' = 6x -> 6x(0) = 0
    std::cout << ca2.derivative().derivative(-10) << std::endl; // x^3'' = 3x^2' = 6x -> 6x(-10) = -60
    std::cout << ca2.derivative().derivative(10) << std::endl; // x^3'' = 3x^2' = 6x -> 6x(0) = 60
    std::cout << std::endl;
    
    std::cout << ca2.integrate(0, 1) << std::endl; // 0.25
    std::cout << ca2.integrate(-1, 1) << std::endl; // 0
    std::cout << ca2.integrate(0, 100) << std::endl; // 25 000 000
    std::cout << ca2.integrate() << std::endl; // 0
    std::cout << std::endl;
    
    try {
        ChebyshevApproximation ca3(sinus, 0, -pi, 1000); // Izuzetak zbog neregularnih granica - Bad parameters
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 0); // Izuzetak zbog neregularnog broja uzoraka - Bad parameters
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.set_m(0); // Izuzetak zbog neregularnog broja uzoraka - Bad order
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.set_m(1001); // Izuzetak zbog neregularnog broja uzoraka - Bad order
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.trunc(-1); // Izuzetak zbog negativnog epsilona - Bad tolerance
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.trunc(1); // Izuzetak zbog prevelikog epsilona (m postaje manje od 1) - Bad tolerance
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3(10); // Izuzetak zbog argumenta izvan granica - Bad argument
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3(-10); // Izuzetak zbog argumenta izvan granica - Bad argument
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.derivative(10); // Izuzetak zbog argumenta izvan granica - Bad argument
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.derivative(-10); // Izuzetak zbog argumenta izvan granica - Bad argument
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.integrate(-10, pi); // Izuzetak zbog granica integracije - Bad interval
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.integrate(pi, -10); // Izuzetak zbog granica integracije - Bad interval
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.integrate(0, 2 * pi); // Izuzetak zbog granica integracije - Bad interval
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        ChebyshevApproximation ca3(sinus, 0, pi, 1000);
        ca3.integrate(2 * pi, 0); // Izuzetak zbog granica integracije - Bad interval
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << ca1.integrate(0, pi/2) << std::endl; // 1
    std::cout << ca1.integrate(pi/2, pi) << std::endl; // 1
    std::cout << ca1.integrate(0, pi) << std::endl; // 2
    std::cout << ca1.integrate() << std::endl; // 2
    std::cout << std::endl;
    
    
    // Testiranje Rombergove metode integracije
    auto par1 = RombergIntegration(sinus, 0, pi/2);
    std::cout << par1.first << " " << par1.second << std::endl; // 1
    auto par2 = RombergIntegration(sinus, pi/2, pi);
    std::cout << par2.first << " " << par2.second << std::endl; // 1
    auto par3 = RombergIntegration(sinus, 0, pi);
    std::cout << par3.first << " " << par3.second << std::endl; // 2
    std::cout << std::endl;
    
    auto par4 = RombergIntegration(kubna, 0, 1);
    std::cout << par4.first << " " << par4.second << std::endl; // 0.25
    auto par5 = RombergIntegration(kubna, -1, 1);
    std::cout << par5.first << " " << par5.second << std::endl; // 0
    auto par6 = RombergIntegration(kubna, 0, 100);
    std::cout << par6.first << " " << par6.second << std::endl; // 25 000 000
    std::cout << std::endl;
    
    try {
        RombergIntegration(kubna, 0, 1, -1); // Izuzetak zbog negativnog epsilona - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        RombergIntegration(kubna, 0, 1, 0.1, -100000, 100); // Izuzetak zbog negativnog nmax - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        RombergIntegration(kubna, 0, 1, 0.1, 100000, -100); // Izuzetak zbog negativnog nmin - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        RombergIntegration(kubna, 0, 1, 0.1, 100, 100000); // Izuzetak zbog nmin > nmax - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    
    // Testiranje integracije metodom smjene tanh sinh
    auto par7 = TanhSinhIntegration(sinus, 0, pi/2);
    std::cout << par7.first << " " << par7.second << std::endl; // 1
    auto par8 = TanhSinhIntegration(sinus, pi/2, pi);
    std::cout << par8.first << " " << par8.second << std::endl; // 1
    auto par9 = TanhSinhIntegration(sinus, 0, pi);
    std::cout << par9.first << " " << par9.second << std::endl; // 2
    std::cout << std::endl;
    
    auto par10 = TanhSinhIntegration(kubna, 0, 1);
    std::cout << par10.first << " " << par10.second << std::endl; // 0.25
    auto par11 = TanhSinhIntegration(kubna, -1, 1);
    std::cout << par11.first << " " << par11.second << std::endl; // 0
    auto par12 = TanhSinhIntegration(kubna, 0, 100);
    std::cout << par12.first << " " << par12.second << std::endl; // 25 000 000
    std::cout << std::endl;
    
    auto fun = [](double x) { return 1 / sqrt(x); };
    auto par13 = TanhSinhIntegration(fun, 0, 1);
    std::cout << par13.first << " " << par13.second << std::endl; // 2
    std::cout << std::endl;
    
    try {
        TanhSinhIntegration(kubna, 0, 1, -1); // Izuzetak zbog negativnog epsilona - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        TanhSinhIntegration(kubna, 0, 1, 0.1, -100000, 100); // Izuzetak zbog negativnog nmax - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        TanhSinhIntegration(kubna, 0, 1, 0.1, 100000, -100); // Izuzetak zbog negativnog nmin - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        TanhSinhIntegration(kubna, 0, 1, 0.1, 100, 100000); // Izuzetak zbog nmin > nmax - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        TanhSinhIntegration(kubna, 0, 1, 0.1, 100000, 100, -5); // Izuzetak zbog negativnog range - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    
    // Testiranje integracije adaptivnom metodom
    auto par15 = AdaptiveIntegration(sinus, 0, pi/2);
    std::cout << par15.first << " " << par15.second << std::endl; // 1
    auto par16 = AdaptiveIntegration(sinus, pi/2, pi);
    std::cout << par16.first << " " << par16.second << std::endl; // 1
    auto par17 = AdaptiveIntegration(sinus, 0, pi);
    std::cout << par17.first << " " << par17.second << std::endl; // 2
    std::cout << std::endl;
    
    auto par18 = AdaptiveIntegration(kubna, 0, 1);
    std::cout << par18.first << " " << par18.second << std::endl; // 0.25
    auto par19 = AdaptiveIntegration(kubna, -1, 1);
    std::cout << par19.first << " " << par19.second << std::endl; // 0
    auto par20 = AdaptiveIntegration(kubna, 0, 100);
    std::cout << par20.first << " " << par20.second << std::endl; // 25 000 000
    std::cout << std::endl;
    
    auto par21 = AdaptiveIntegration(fun, 0, 1);
    std::cout << par21.first << " " << par21.second << std::endl; // 2
    std::cout << std::endl;
    
    try {
        TanhSinhIntegration(kubna, 0, 1, -1); // Izuzetak zbog negativnog epsilona - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        TanhSinhIntegration(kubna, 0, 1, 0.1, -5); // Izuzetak zbog negativne max dubine - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        TanhSinhIntegration(kubna, 0, 1, 0.1, 50, -5); // Izuzetak zbog negativnog nmin - Bad parameter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
}
