#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>
#include <stdexcept>

int max(int a, int b) {
    if(a > b) return a;
    else return b;
}

int min(int a, int b) {
    if(a < b) return a;
    else return b;
}

class AbstractInterpolator {
    protected:
    std::vector<std::pair<double, double>> parovi;
    double eps = std::numeric_limits<double>::epsilon();
    mutable int interval;
    int Locate(double x) const;
    public:
    AbstractInterpolator(const std::vector<std::pair<double, double>> &data);
    virtual double operator()(double x) const = 0;
};

AbstractInterpolator::AbstractInterpolator(const std::vector<std::pair<double, double>> &data) {
    parovi = data;
    interval = 0;
    std::sort(parovi.begin(), parovi.end(), [](std::pair<double, double> prvi, std::pair<double, double> drugi) {
        return prvi.first < drugi.first;
    });
    for(unsigned int i(0); i < parovi.size() - 1; i++)
        if(std::abs(parovi[i].first - parovi[i + 1].first) < eps) throw std::domain_error("Invalid data set");
}

int AbstractInterpolator::Locate(double x) const {
    if(x < parovi[0].first || std::abs(std::abs(x) - std::abs(parovi[0].first)) < eps) return 0;
    if(x > parovi[parovi.size() - 1].first) return parovi.size();
    int dno(0), vrh(parovi.size() - 1), srednji;
    while(vrh >= dno) {
        srednji = (dno + vrh) / 2;
        if(std::abs(std::abs(parovi[srednji].first) - std::abs(x)) < eps) return srednji;
        else if(parovi[srednji].first > x) vrh = srednji - 1;
        else dno = srednji + 1;
    }
    if(parovi[srednji].first > x) return srednji;
    else return srednji + 1;
}

class LinearInterpolator : public AbstractInterpolator {
    public:
    LinearInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {}
    double operator()(double x) const override;
};

double LinearInterpolator::operator()(double x) const {
    interval = Locate(x);
    if(interval == 0)
        return ((parovi[1].first - x) / (parovi[1].first - parovi[0].first)) * parovi[0].second +
                (x - (parovi[0].first) / (parovi[1].first - parovi[0].first)) * parovi[1].second;
    else if(interval == parovi.size()) {
        int novi(parovi.size() - 2);
        return parovi[novi].second + ((parovi[novi + 1].second - parovi[novi].second) / (parovi[novi + 1].first - parovi[novi].first)) * (x - parovi[novi].first);
        
    }
    else {
        int novi(interval - 1);
        return ((parovi[novi + 1].first - x) / (parovi[novi + 1].first - parovi[novi].first)) * parovi[novi].second +
                (x - (parovi[novi].first) / (parovi[novi + 1].first - parovi[novi].first)) * parovi[novi + 1].second;

    }        
}

class PolynomialInterpolator : public AbstractInterpolator {
    std::vector<double> koeficijenti;
    std::vector<double> koeficijenti2;
    public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data);
    double operator()(double x) const override;
    void AddPoint(const std::pair<double, double> &p);
    std::vector<double> GetCoefficients() const;
    
    void Ispisii() {
        for(int i(0); i < koeficijenti.size(); i++) {
            std::cout << koeficijenti[i] << " ";
        }
    }
};

PolynomialInterpolator::PolynomialInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {
    int n(parovi.size() - 1);
    koeficijenti.resize(parovi.size());
    for(int i(0); i < koeficijenti.size(); i++)
        koeficijenti[i] = parovi[i].second;
    
    koeficijenti2.push_back(parovi[n].first);
    koeficijenti2.push_back(parovi[n].second);
    
    for(int j(1); j <= n; j++)
        for(int i(n); i >= j; i--) {
            koeficijenti[i] = (koeficijenti[i] - koeficijenti[i - 1]) / (parovi[i].first - parovi[i - j].first);
            if(i == n) koeficijenti2.push_back(koeficijenti[i]);
        }
}

inline double PolynomialInterpolator::operator()(double x) const {
    int n(koeficijenti.size() - 1);
    double f(koeficijenti[n]);
    for(int i(n); i >= 1; i--)
        f = f * (x - parovi[i - 1].first) + koeficijenti[i - 1];
    return f;
}

void PolynomialInterpolator::AddPoint(const std::pair<double, double> &p) {
    for(int i(0); i < parovi.size(); i++)
        if(std::abs(p.first - parovi[i].first) < eps) throw std::domain_error("Invalid point");
    parovi.push_back(p);
    std::vector<double> pom(koeficijenti2.size() + 1);
    pom[0] = p.first;
    pom[1] = p.second;
    pom[2] = (pom[1] - parovi[parovi.size() - 2].second) / (pom[0] - parovi[parovi.size() - 2].first); 
    for(int i(2); i < pom.size(); i++)
    {
        pom[i] = (pom[i - 1] - koeficijenti2[i - 1]) / (parovi[parovi.size() - 1].first - parovi[parovi.size() - i].first);
    }
    koeficijenti.push_back(pom[pom.size() - 1]);
    koeficijenti2 = pom;
}

std::vector<double> PolynomialInterpolator::GetCoefficients() const {
    std::vector<double> w(parovi.size() + 1);
    std::vector<double> p(parovi.size());
    int n(parovi.size() - 1);
    w[0] = 1;
    for(int i(1); i <= n + 1 ; i++) {
        w[i] = w[i - 1];
        for(int j(i - 1); j >= 1; j--)
            w[j] = w[j - 1] - parovi[i - 1].first * w[j];
        w[0] = -parovi[i - 1].first * w[0];
    }
    for(int i(0); i <= n; i++) {
        double f = 1;
        for(int j(0); j <= n; j++)
            if(j != i)
                f = f * (parovi[i].first - parovi[j].first);
        f = parovi[i].second / f;
        std::vector<double> v(parovi.size() + 1);
        for(int j(0); j <= n + 1; j++)
            v[j] = w[j];
        for(int j(n); j >= 0; j--) {
            v[j] = v[j] + parovi[i].first * v[j + 1];
            p[j] = p[j] + f * v[j + 1];
        }
    }
    return p;
}


class PiecewisePolynomialInterpolator : public AbstractInterpolator {
    int k;
    public:
    PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data, int order) : AbstractInterpolator(data) {
        if(order < 1 || order > data.size()) throw std::domain_error("Invalid order");
        k = order;
    }
    double operator()(double x) const override;
};

double PiecewisePolynomialInterpolator::operator()(double x) const {
    interval = Locate(x);
    int j1, j2;
    if(k % 2 == 0) {
        j1 = interval - 1 - k / 2;
        j2 = interval - 1 + k / 2;
        //td::cout << "-1-";
        if(j1 < 0) {
            //std::cout << "-2-";
            j1 = 0;
            j2 = k + 1;
        }
        if(j2 > parovi.size()) {
            //std::cout << "-3-";
            j1 = parovi.size() - 1 - k;
            j2 = parovi.size() - 1;
        }
    }
    else {
        j1 = interval - 1 - (k - 1) / 2;
        j2 = interval - 1 + (k + 1) / 2;
        //std::cout << "-4-"; 
        if(j1 < 0) {
            //std::cout << "-5-";
            j1 = 0;
            j2 = k;
        }
        if(j2 >= parovi.size()) {
            //std::cout << "-6-";
            j1 = parovi.size() - 1 - k;
            j2 = parovi.size() - 1;
        }
    }
    
    double s(0);
    for(int i(j1); i <= j2; i++) {
        double p(parovi[i].second);
        for(int j(j1); j <= j2; j++)
            if(j != i)
                p = p * (x - parovi[j].first) / (parovi[i].first - parovi[j].first);
        s = s + p;
    }
    return s;
}


class BarycentricInterpolator : public AbstractInterpolator {
    int k;
    std::vector<double> koeficijenti;
    public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order);
    double operator()(double x) const override;
    std::vector<double> GetWeights() const { return koeficijenti; }
};

BarycentricInterpolator::BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order)
    : AbstractInterpolator(data) {
    if(order < 0 || order > data.size()) throw std::domain_error("Invalid order");
    
        int n(parovi.size());
        koeficijenti.resize(n);
        
        for(int i(0); i < n; i++) {
            koeficijenti[i] = 0;
            double p;
            for(int k = max(1, i - order); k >= min(i, n - order); k--) {
                p = 1;
                for(int j(k); j < k + order; j++)
                    if(j != i)
                        p = p / (parovi[i].first - parovi[j].first);
                if(k % 2 != 0) 
                    p = -p;
            }
            koeficijenti[i] += p;
        }
        
}

double BarycentricInterpolator::operator()(double x) const {
    double p(0);
    double q(0);
    for(int i(0); i < parovi.size(); i++) {
        if(std::abs(x - parovi[i].first) < eps)
            return parovi[i].second;
        double u(koeficijenti[i] / (x - parovi[i].first));
        p = p + u * parovi[i].second;
        q = q + u;
    }
    return p / q;
}


template <typename FunTip>
double Limit(FunTip f, double x0, double h = 0, double eps = 1e-8, double nmax = 20) {
    if(eps <= 0) throw std::domain_error("Invalid parameters");
    if(nmax < 3 || nmax > 30) throw std::domain_error("Invalid parameters");
    if(std::abs(h) < eps)
        h = 0.001 * max(1, std::abs(x0));
        
    double old = std::numeric_limits<double>::infinity();
    std::vector<double> y(nmax + 1);
    double p;
    int i(1);
    for(;;) {
        if(i > nmax) throw std::logic_error("Accuracy goal is not achieved");
        y[i] = f(x0 + h);
        p = 2;
        for(int k(i - 1); k >= 1; k--) {
            y[k] = (p * y[k + 1] - y[k]) / (p - 1);
            p = 2 * p;
        }
        if(std::abs(y[1] - old) < eps) return y[1];
        old = y[1];
        h = h / 2;
        i++;
    }
    return y[1];
}


int main ()
{
    // Testiranje konstruktora bazne funkcije
    /*
    std::vector<std::pair<double, double>> v;
    for(double i(0); i < 10; i += 0.4) {
        std::pair<double, double> p;
        p.first = i * 10.8 / 4;
        p.second = i;
        v.push_back(p);
    }
    AbstractInterpolator ai(v);
    ai.Ispisi();
    */
    
    // Testiranje funkcije Locate
    /*
    std::vector<std::pair<double, double>> v1;
    for(double i(0); i < 10; i++) {
        std::pair<double, double> p;
        p.first = i * 2;
        p.second = i;
        v1.push_back(p);
    }
    AbstractInterpolator ai1(v1);
    ai1.Ispisi();
    std::cout << ai1.Locate(5) << std::endl;
    std::cout << ai1.Locate(0) << std::endl;
    std::cout << ai1.Locate(20) << std::endl;
    */
    
    // Testiranje linearne interpolacije
    std::vector<std::pair<double, double>> v2;
    // Pravac
    for(double i(0); i < 10; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i;
        v2.push_back(p);
    }
    LinearInterpolator li(v2);
    std::cout << li(1) << std::endl;
    std::cout << li(5) << std::endl;
    std::cout << li(6) << std::endl;
    std::cout << li(-5) << std::endl;
    std::cout << li(10) << std::endl;
    std::cout << li(100) << std::endl;
    
    v2[5].first = 0;
    v2[5].second = 0;
    try {
        LinearInterpolator li1(v2);
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    // x kvadrat
    v2.clear();
    for(double i(-10000); i < 10000; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i * i;
        v2.push_back(p);
    }
    LinearInterpolator li1(v2);
    std::cout << li1(1) << std::endl;
    std::cout << li1(5) << std::endl;
    std::cout << li1(6) << std::endl;
    std::cout << li1(-5) << std::endl;
    std::cout << li1(10) << std::endl;
    std::cout << li1(100) << std::endl;
    
    
    // Testiranje polinomske interpolacije
    std::vector<std::pair<double, double>> v3;
    // Pravac
    for(double i(0); i < 10; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i;
        v3.push_back(p);
    }
    PolynomialInterpolator pi(v3);
    std::cout << std::endl;
    std::cout << pi(1) << std::endl;
    std::cout << pi(5) << std::endl;
    std::cout << pi(6) << std::endl;
    std::cout << pi(-5) << std::endl;
    std::cout << pi(10) << std::endl;
    std::cout << pi(100) << std::endl;
    
    v3[5].first = 0;
    v3[5].second = 0;
    try {
        PolynomialInterpolator pi1(v3);
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    pi.AddPoint({11, 11});
    try {
        pi.AddPoint({11, 11});
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    // x kvadrat
    v3.clear();
    for(double i(-1000); i < 1000; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i * i;
        v3.push_back(p);
    }
    PolynomialInterpolator pi1(v3);
    std::cout << std::endl;
    std::cout << pi1(1) << std::endl;
    std::cout << pi1(5) << std::endl;
    std::cout << pi1(6) << std::endl;
    std::cout << pi1(-5) << std::endl;
    std::cout << pi1(10) << std::endl;
    std::cout << pi1(100) << std::endl;
    
    // Vektor
    std::vector<double> p(pi.GetCoefficients());
    for(unsigned int i(0); i < p.size(); i++)
        std::cout << p[i] << " ";
    
    
    // Testiranje interpolacije polinomima dio po dio
    std::vector<std::pair<double, double>> v4;
    // Pravac
    for(double i(0); i < 10; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i;
        v4.push_back(p);
    }
    PiecewisePolynomialInterpolator pw(v4, 4);
    std::cout << std::endl;
    std::cout << pw(1) << std::endl;
    std::cout << pw(5) << std::endl;
    std::cout << pw(6) << std::endl;
    std::cout << pw(-5) << std::endl;
    std::cout << pw(10) << std::endl;
    
    v4[5].first = 0;
    v4[5].second = 0;
    try {
        PiecewisePolynomialInterpolator pw1(v4, 4);
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    // x kvadrat
    v4.clear();
    for(double i(-1000); i < 1000; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i * i;
        v4.push_back(p);
    }
    PiecewisePolynomialInterpolator pw1(v4, 3);
    std::cout << std::endl;
    std::cout << pw1(1) << std::endl;
    std::cout << pw1(5) << std::endl;
    std::cout << pw1(6) << std::endl;
    std::cout << pw1(-5) << std::endl;
    std::cout << pw1(10) << std::endl;
    std::cout << pw1(100) << std::endl;
    
    
    // Testiranje baricentricne interpolacije
    std::vector<std::pair<double, double>> v5;
    // Pravac
    for(double i(0); i < 10; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i;
        v5.push_back(p);
    }
    BarycentricInterpolator bi(v5, 3);
    std::cout << std::endl;
    std::cout << bi(1) << std::endl;
    std::cout << bi(5) << std::endl;
    std::cout << bi(6) << std::endl;
    std::cout << bi(-5) << std::endl;
    std::cout << bi(10) << std::endl;
    
    v5[5].first = 0;
    v5[5].second = 0;
    try {
        BarycentricInterpolator bi1(v5, 3);
    }
    catch(std::domain_error izuzetak) {
        std::cout << izuzetak.what() << std::endl;
    }
    
    // x kvadrat
    v5.clear();
    for(double i(-1000); i < 1000; i++) {
        std::pair<double, double> p;
        p.first = i;
        p.second = i * i;
        v5.push_back(p);
    }
    PiecewisePolynomialInterpolator bi1(v5, 3);
    std::cout << std::endl;
    std::cout << bi1(1) << std::endl;
    std::cout << bi1(5) << std::endl;
    std::cout << bi1(6) << std::endl;
    std::cout << bi1(-5) << std::endl;
    std::cout << bi1(10) << std::endl;
    std::cout << bi1(100) << std::endl;

    //std::cout<<Limit([](double x) {return std::pow(x,1./3);},0,0,0,5)<<std::endl;
    std::cout<<Limit([](double x) {return std::log(x)/(x-1);},1);
    std::cout<<Limit([](double x) {return 1 / x;},1);

   return 0;
}