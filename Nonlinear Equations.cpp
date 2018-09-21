#include <iostream>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <complex>
#include <vector>
#include <utility>

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template <typename FunTip>
bool BracketRoot(FunTip f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
    if(hinit < 0 || hmax < 0 || lambda < 0) throw std::domain_error("Invalid parameters");
    int brojac(0);
    double hinitcopy(hinit);
    
    pocetak:
        a = x0;
        double f1(f(a));
        
        while(std::abs(hinit) < hmax) {
            b = a + hinit;
            double f2(f(b));
            if(f1 * f2 <= 0) {
                if (b < a) {
                    double temp(a);
                    a = b; 
                    b = temp;
                }
                return true;
            }
            hinit *= lambda;
            a = b;
            f1 = f2;
        }
    
    if(brojac == 0) {
        x0 *= -1;
        hinit = hinitcopy;
        brojac++;
        goto pocetak;
    }
    
    return false;
}

template <typename FunTip>
double RegulaFalsiSolve(FunTip f, double a, double b, RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100) {
    if((f(a) > 0 && f(b) > 0) || (f(a) < 0 && f(b) < 0)) throw std::range_error("Root must be bracketed");
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    
    int iteracije(0);
    if(mode == Unmodified) {
        double f1(f(a));
        double f2(f(b));
        double c(a);
        double cold(b);
        
        while(std::abs(c - cold) > eps) {
            if(iteracije > maxiter) throw std::logic_error("Given accuracy has not achieved");
            cold = c;
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3(f(c));
            if(f3 == 0) return c;
            if(f1 * f3 < 0) {
                b = c;
                f2 = f3;
            }
            else {
                a = c;
                f1 = f3;
            }
            iteracije++;
        }
        
        return c;
    }
    else if(mode == Illinois) {
        double f1(f(a));
        double f2(f(b));
        double c(a);
        
        while(std::abs(b - a) > eps) {
            if(iteracije > maxiter) throw std::logic_error("Given accuracy has not achieved");
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3(f(c));
            if(f3 == 0) return c;
            if(f1 * f3 < 0) {
                b = a;
                f2 = f1;
            }
            else 
                f2 /= 2;
            a = c;
            f1 = f3;
            iteracije++;
        }
        
        return c;
    }
    else if(mode == Slavic) {
        auto phi = [f](double x) { return f(x) / (1 + std::abs(f(x))); }; // Izmijenjena funkcija f
        double f1(phi(a));
        double f2(phi(b));
        double c(a);
        double cold(b);
        
        while(std::abs(c - cold) > eps) {
            if(iteracije > maxiter) throw std::logic_error("Given accuracy has not achieved");
            cold = c;
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3(f(c));
            if(f3 == 0) return c;
            if(f1 * f3 < 0) {
                b = c;
                f2 = f3;
            }
            else {
                a = c;
                f1 = f3;
            }
            iteracije++;
        }
        
        return c;
    }
    else if(mode == IllinoisSlavic) {
        auto phi = [f](double x) { return f(x) / (1 + std::abs(f(x))); }; // Izmijenjena funkcija f
        double f1(phi(a));
        double f2(phi(b));
        double c(a);
        
        while(std::abs(b - a) > eps) {
            if(iteracije > maxiter) throw std::logic_error("Given accuracy has not achieved");
            c = (a * f2 - b * f1) / (f2 - f1);
            double f3(f(c));
            if(f3 == 0) return c;
            if(f1 * f3 < 0) {
                b = a;
                f2 = f1;
            }
            else 
                f2 /= 2;
            a = c;
            f1 = f3;
            iteracije++;
        }
        
        return c;
    }
    else {}
}

template <typename FunTip>
double RiddersSolve(FunTip f, double a, double b, double eps = 1e-10, int maxiter = 100) {
    if((f(a) > 0 && f(b) > 0) || (f(a) < 0 && f(b) < 0)) throw std::range_error("Root must be bracketed");
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    
    auto sgn = [](double x) {
        if(x < 0) return - 1;
        else if(x > 0) return 1;
        else return 0;
    };
    int iteracije(0);
    double f1(f(a));
    double f2(f(b));
    
    while(std::abs(b - a) > eps) {
        if(iteracije > maxiter) throw std::logic_error("Given accuracy has not achieved");
        double c((a + b) / 2);
        double f3(f(c));
        if(f3 == 0) return c;
        double d(c + f3 * (c - a) * sgn(f1 - f2) / std::sqrt(f3 * f3 - f1 * f2));
        double f4 = f(d);
        if(f4 == 0) return d;
        if(f3 * f4 <= 0) {
            a = c;
            b = d;
            f1 = f3;
            f2 = f4;
        }
        else if(f1 * f4 <= 0) {
            b = d;
            f2 = f4;
        }
        else {
            a = d;
            f1 = f4;
        }
        iteracije++;
    }
    
    return (a + b) / 2;
}

template <typename FunTip1, typename FunTip2>
double NewtonRaphsonSolve(FunTip1 f, FunTip2 fprim, double x0, double eps = 1e-10, int maxiter = 100) {
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    
    int iteracije(0);
    double deltax(std::numeric_limits<double>::infinity());
    while(std::abs(deltax) > eps) {
        if(iteracije > maxiter) throw std::logic_error("Convergence has not achieved");
        double v(f(x0));
        if(v == 0) return x0;
        if(fprim(x0) == 0) throw std::logic_error("Convergence has not achieved");
        deltax = v / fprim(x0);
        x0 -= deltax;
        iteracije++;
    }
    
    return x0;
}

std::pair<std::complex<double>, bool> Laguerre(std::vector<double> coefficients, int n, std::complex<double> x, double eps = 1e-10,
int maxtrials = 10) {
    std::complex<double> deltax(std::numeric_limits<double>::infinity());
    int k(1);
    while(std::abs(deltax) > eps && k < maxtrials) {
        std::complex<double> f(coefficients[n]);
        std::complex<double> d;
        d.real(0);
        d.imag(0);
        std::complex<double> s;
        s.real(0);
        s.imag(0);
        for(int i(n - 1); i >= 0; i--) {
            s = s * x + std::complex<double>(2) * d;
            d = d * x + f;
            f = f * x + coefficients[i];
        }
        if(f == std::complex<double>(0)) return std::pair<std::complex<double>, bool>(x, true);
        std::complex<double> r(std::sqrt(std::complex<double>(n - 1) * (std::complex<double>(n - 1) * d * d - std::complex<double>(n) * f * s)));
        if(std::abs(d + r) > std::abs(d - r))
            deltax = std::complex<double>(n) * f / (d + r);
        else
            deltax = std::complex<double>(n) * f / (d - r);
        x -= deltax;
        k++;
    }
    
    if(std::abs(deltax) <= eps) return std::pair<std::complex<double>, bool>(x, true);
    return std::pair<std::complex<double>, bool>(x, false);
}

std::vector<std::complex<double>>PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10,
int maxiters = 100, int maxtrials = 10) {
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    
    for(int i(coefficients.size() - 1); i >= 1; i--) {
        int t(1);
        bool c(false);
        while(!c && t < maxtrials) {
            int re((rand() % 20) - 10), im((rand() % 20) - 10);
            std::complex<double> x;
            x.real(re);
            x.imag(im);
            
        }
    }
}

std::vector< std::complex<double>> PolyRoots(std::vector<double> coefficients, double eps = 1e-10,
int maxiters = 100, int maxtrials = 10) {
    int n(coefficients.size() - 1);
    int i(n);
    std::vector< std::complex<double>> rez(n);
    
    while(i >= 1) {
        int t(1);
        bool c(false);
        std::complex<double> x;
        while(!c && t < maxtrials) {
            int re((rand() % 20) - 10), im((rand() % 20) - 10);
            x.real(re);
            x.imag(im);
            std::pair<std::complex<double>, bool> par(Laguerre(coefficients, i, x));
            x = par.first;
            c = par.second;
            t++;
        }
        if(!c) throw std::logic_error("Convergence has not achieved");
        if(std::abs(std::imag(x)) < eps) {
            double xx = std::real(x);
            rez[i] = x;
            double v(coefficients[i]);
            for(int j(i - 1); j >= 0; j--) {
                double w(coefficients[j]);
                coefficients[j] = v;
                v = w + xx * v;
            }
            i--;
        }
        else {
            rez[i] = std::complex<double>(x);
            rez[i - 1] = std::conj(x);
            double alfa(2 * std::real(x));
            double beta(std::abs(x) * std::abs(x));
            double u(coefficients[i]);
            double v(coefficients[i - 1] + alfa * u);
            for(int j(i - 1); j >= 0; j--) {
                double w(coefficients[j]);
                coefficients[j] = u;
                u = v;
                v = w + alfa * v - beta * coefficients[j];
            }
            i -= 2;
        }
    }
    return rez;
}

int main () {
    
    // Testiranje BracketRoot funkcije
    auto f1 = [](double x) { return x * x; }; // kvadratna - funkcija samo dira x - osu tako da će biti false
    double a1, b1;
    if(BracketRoot(f1, 1, a1, b1)) std::cout << a1 << " " << b1 << std::endl;
    else std::cout << "Nije pronadjena nula!" << std::endl;
    
    auto f2 = [](double x) { return log(x); }; // ln(x) - nula je u x = 1
    double a2, b2;
    if(BracketRoot(f2, 1, a2, b2)) std::cout << a2 << " " << b2 << std::endl;
    else std::cout << "Nije pronadjena nula!" << std::endl;
    
    auto f3 = [](double x) { return x; }; // x - nula je u x = 0
    double a3, b3;
    if(BracketRoot(f3, 1, a3, b3)) std::cout << a3 << " " << b3 << std::endl;
    else std::cout << "Nije pronadjena nula!" << std::endl;
    
    auto f4 = [](double x) { return (x * x - 10) / 2; }; // nula je u x = 3.16228
    double a4, b4;
    if(BracketRoot(f4, 1, a4, b4)) std::cout << a4 << " " << b4 << std::endl;
    else std::cout << "Nije pronadjena nula!" << std::endl;
    
    auto f5 = [](double x) { return std::sin(x); }; // nula je u 0 ili pi, zavisi koja se početna tačka proslijedi
    double a5, b5;
    if(BracketRoot(f5, -1, a5, b5)) std::cout << a5 << " " << b5 << std::endl;
    else std::cout << "Nije pronadjena nula!" << std::endl;
    
    // Izuzeci
    try {
        BracketRoot(f1, 1, a1, b1, -1); // Negativan hinit
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        BracketRoot(f1, 1, a1, b1, 1e-5, -1); // Negativan hmax
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        BracketRoot(f1, 1, a1, b1, 1e-5, 1e10, -1); // Negativan lambda
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    
    
    // Testiranje RegulaFalsi funkcije
    // Unmodified
    try {
        std::cout << RegulaFalsiSolve(f2, a2, b2, RegulaFalsiMode::Unmodified) << std::endl; // 1
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
    std::cout << RegulaFalsiSolve(f3, a3, b3, RegulaFalsiMode::Unmodified) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f4, a4, b4, RegulaFalsiMode::Unmodified) << std::endl; // 3.16228
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f5, a5, b5, RegulaFalsiMode::Unmodified) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    // Illinois
    try {
        std::cout << RegulaFalsiSolve(f2, a2, b2, RegulaFalsiMode::Illinois) << std::endl; // 1
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f3, a3, b3, RegulaFalsiMode::Illinois) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f4, a4, b4, RegulaFalsiMode::Illinois) << std::endl; // 3.16228
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f5, a5, b5, RegulaFalsiMode::Illinois) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    // Slavic
    try {
        std::cout << RegulaFalsiSolve(f2, a2, b2, RegulaFalsiMode::Slavic) << std::endl; // 1
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f3, a3, b3, RegulaFalsiMode::Slavic) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f4, a4, b4, RegulaFalsiMode::Slavic) << std::endl; // 3.16228
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f5, a5, b5, RegulaFalsiMode::Slavic) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    // IllinoisSlavic
    try {
        std::cout << RegulaFalsiSolve(f2, a2, b2, RegulaFalsiMode::IllinoisSlavic) << std::endl; // 1
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f3, a3, b3, RegulaFalsiMode::IllinoisSlavic) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f4, a4, b4, RegulaFalsiMode::IllinoisSlavic) << std::endl; // 3.16228
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        std::cout << RegulaFalsiSolve(f5, a5, b5, RegulaFalsiMode::IllinoisSlavic) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    // Izuzeci
    try {
        RegulaFalsiSolve(f1, a2, b2); // f(a) i f(b) istog znaka
    }
    catch(std::range_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        RegulaFalsiSolve(f2, a2, b2, RegulaFalsiMode::Slavic, -1, 100); // Negativan eps
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        RegulaFalsiSolve(f2, a2, b2, RegulaFalsiMode::Slavic, 1, -100); // Negativan maxiter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    std::cout << std::endl;
    
    
    
    // Testiranje NewtonRaphsonSolve funkcije
    auto fprim2 = [](double x) { return 1 / x; };
    auto fprim3 = [](double x) { return 1; };
    auto fprim4 = [](double x) { return x; };
    auto fprim5 = [](double x) { return std::cos(x); };
    
    try {
        std::cout << NewtonRaphsonSolve(f2, fprim2, 1) << std::endl; // 1
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
    std::cout << NewtonRaphsonSolve(f3, fprim3, 1) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
    std::cout << NewtonRaphsonSolve(f4, fprim4, 1) << std::endl; // 3.16228
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
    std::cout << NewtonRaphsonSolve(f5, fprim5, 1) << std::endl; // 0
    }
    catch(std::logic_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    // Izuzeci
    try {
        NewtonRaphsonSolve(f2, fprim2, 1, -1); // Negativan eps
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    try {
        NewtonRaphsonSolve(f2, fprim2, 1, 1, -1); // Negativan maxiter
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    


    
	return 0;
}