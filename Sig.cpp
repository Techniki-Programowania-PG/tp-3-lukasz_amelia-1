#include <cmath>
#include <matplot/matplot.h>
#include <vector>
#include <complex>
#include <iostream>
#include <numeric>
#include <pybind11/pybind11.h>


using namespace matplot;
namespace py = pybind11;

void plotSin(double amplitude, int frequency, double begin, double end, int samplesCount) {
    std::vector<double> x = linspace(begin, end, samplesCount);

    std::vector<double> y = transform(x, [&](auto x) {
        double k = frequency * x;
        return (amplitude * sin(k)); });
    plot(x, y);

    show();
}

void plotCos(double amplitude, int frequency, double begin, double end, int samplesCount) {
    std::vector<double> x = linspace(begin, end, samplesCount);

    std::vector<double> y = transform(x, [&](auto x) {
        double k = frequency * x;
        return (amplitude * cos(k)); 
    });
    plot(x, y);

    show();
}

void plotSquare(double amplitude, int frequency, double begin, double end, int samplesCount) {
   std::vector<double> x = linspace(begin, end, samplesCount);

   std::vector<double> y = transform(x, [&](auto x) {
       double k = frequency * x;
       if (sin(k) > 0) {
           return amplitude;
       }
       else {
           return 0.0;
       }
            
            });
    plot(x, y);

    show();
    
}

void plotSaw(double amplitude, int frequency, double begin, double end) {
    std::vector<double> x;

    for (int i = 0; i < end * frequency; i += frequency) {
        x.push_back(i); x.push_back(i);
    }

    std::vector<double> y1;
    for (int i = 0; i < end; i++) {
        y1.push_back(0);
        y1.push_back(amplitude);
    }
    
    plot(x, y1);

    show();

}

std::vector<std::complex<double>> DFT(const std::vector<double>& x) {
    size_t N = x.size();
    std::vector<std::complex<double>> X(N);

    for (size_t k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t n = 0; n < N; ++n) {
            double angle = -2 * pi * k * n / N;
            sum += x[n] * std::complex<double>(cos(angle), sin(angle));
        }
        X[k] = sum;
    }
    return X;
}

std::vector<double> IDFT(const std::vector<std::complex<double>>& X) {
    size_t N = X.size();
    std::vector<std::complex<double>> x_complex(N);

    for (size_t n = 0; n < N; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t k = 0; k < N; ++k) {
            double angle = 2 * pi * k * n / N;
            sum += X[k] * std::complex<double>(cos(angle), sin(angle));
        }
        x_complex[n] = sum / static_cast<double>(N);
    }


    std::vector<double> x_real(N);
    for (size_t i = 0; i < N; ++i) {
        x_real[i] = x_complex[i].real();
    }

    return x_real;
}

std::vector<double> filterOneD(const std::vector<double>& signal, size_t windowSize) {

    std::vector<double> filtered;
    size_t n = signal.size();
    for (size_t i = 0; i <= n - windowSize; ++i) {
        double sum = accumulate(signal.begin() + i, signal.begin() + i + windowSize, 0.0);
        filtered.push_back(sum / windowSize);
    }

    return filtered;
}

std::vector<std::vector<double>> filterTwoD(const std::vector<std::vector<double>>& input) {
    size_t rows = input.size();
    size_t cols = input[0].size();
    std::vector<std::vector<double>> output(rows - 2, std::vector<double>(cols - 2, 0.0));

    for (size_t i = 1; i < rows - 1; ++i) {
        for (size_t j = 1; j < cols - 1; ++j) {
            double sum = 0.0;
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    sum += input[i + dx][j + dy];
                }
            }
            output[i - 1][j - 1] = sum / 9.0;
        }
    }
    return output;
}

std::vector<double> derivative(const std::vector<double>&y, double dx) {
    std::vector<double> dy;
    for (size_t i = 0; i < y.size() - 1; ++i) {
        dy.push_back((y[i + 1] - y[i]) / dx);
    }
    return dy;
}


void plotDFT() {
    std::cout << "Dokonuje Dyskretnej Transformaty Fouriera na przykładowym sygnale - 64 probki, 30 Hz. Wyswietlony zostanie najpierw sygnal oryginalny, a potem jego Dyskretna Transformata Fouriera" << std::endl;
    int N = 64;
    std::vector<double> t(N), signal(N);
    for (int n = 0; n < N; ++n) {
        t[n] = n;

        signal[n] = sin(2 * pi * n / N * 30);
    }


    std::vector<std::complex<double>> X = DFT(signal);
    std::cout << "Dokonuje Odwrotnej Dyskretnej Transformaty Fouriera" << std::endl;
    std::vector<double> magnitude(N);
    for (int i = 0; i < N; ++i) {
        magnitude[i] = abs(X[i]);
    }

    plot(t, signal);
    show();

    stem(t, magnitude);
    show();

    std::vector<double> recX = IDFT(X);
    plot(t, recX);
    show();
}

void plotFiltering() {
    std::cout << "Sygnal oryginalny: " << std::endl;
    int amplitude = 2;
    int frequency = 5;
    int begin = 0;
    int end = 2;
    int samplesCount = 100;
    std::vector<double> x = linspace(begin, end, samplesCount);
    int i = 0;

    std::vector<double> y = transform(x, [&](auto x) {
        i++;
        double k = frequency * x;
        if (i % 4 == 0) {
            return (amplitude * sin(k) + 0.05);
        }
        else {
            return (amplitude * sin(k) + 0.1);
        }
        });
    plot(x, y);

    show();

    size_t window_size = 5;
    std::vector<double> y_filtered = filterOneD(y, window_size);
    std::vector<double> x_filtered(x.begin() + window_size - 1, x.end());

    std::cout << "Sygnal przefiltrowany: " << std::endl;

    plot(x_filtered, y_filtered);
    show();

    int rows = 75;
    int cols = 75;
    std::vector<std::vector<double>> image(rows, std::vector<double>(cols));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (i % 2 == 0) {
                double value = sin(i * 0.1) + cos(j * 0.15) + 0.1;
                image[i][j] = value;
            }
            else {
                double value = sin(i * 0.1) + cos(j * 0.15) + 0.2;
                image[i][j] = value;
            }
        }
    }

    std::cout << "Oryginalny obraz" << std::endl;
    imagesc(image);
    show();

    std::vector<std::vector<double>> filtered = filterTwoD(image);

    std::cout << "Przefiltrowany obraz" << std::endl;
    imagesc(filtered);
    show();
}

void plotDerivative() {
    std::cout << "Sygnal oryginalny: sin(x)" << std::endl;
    int n = 100;
    double sample = 0.1;

    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i * sample;
        y[i] = sin(x[i]);
    }

    plot(x, y);
    show();

    std::vector<double> yprim = derivative(y, sample);
    std::vector<double> yprimx(x.begin(), x.end());
    std::cout << "Pochodna sygnalu: " << std::endl;
    plot(yprimx, yprim);
    show();
}

PYBIND11_MODULE(Sig, m) {
    m.doc() = "test";
    m.def("plotSin", &plotSin, "sinus",
        pybind11::arg("amplitude"), pybind11::arg("frequency"), pybind11::arg("begin"), pybind11::arg("end"), pybind11::arg("samplesCount"));
    m.def("plotCos", &plotCos, "cosinus",
        pybind11::arg("amplitude"), pybind11::arg("frequency"), pybind11::arg("begin"), pybind11::arg("end"), pybind11::arg("samplesCount"));
    m.def("plotSquare", &plotSquare, "kwadrat",
        pybind11::arg("amplitude"), pybind11::arg("frequency"), pybind11::arg("begin"), pybind11::arg("end"), pybind11::arg("samplesCount"));
    m.def("plotSaw", &plotSaw, "pila",
        pybind11::arg("amplitude"), pybind11::arg("frequency"), pybind11::arg("begin"), pybind11::arg("end"));
    m.def("DFT", &DFT, "dyskretna transformata fouriera",
        pybind11::arg("x"));
    m.def("IDFT", &IDFT, "odwrotna dyskretna transformata fouriera",
        pybind11::arg("x"));
    m.def("filterOneD", &filterOneD, "filtracja 1D",
        pybind11::arg("signal"), pybind11::arg("windowSize"));
    m.def("filterTwoD", &filterTwoD, "filtracja 2D",
        pybind11::arg("input"));
    m.def("plotDFT", &plotDFT, "transformata fouriera");
    m.def("plotFiltering", &plotFiltering, "filtrowanie");
    m.def("plotDerivative", &plotDerivative, "pochodna sygnalu");
}

int main() {

    std::cout << "Program do analizy sygnalow automatyki. Dostepne opcje: \n 1 - Dyskretna Transformata Fouriera - normalna i odwrotna \n 2 - Filtracja 1D \n 3 - Generator sygnalow \n 4 - Pochodna sygnalu" << std::endl;
    int choice;
    std::cin >> choice;
    switch (choice) {
    case 1: {
        

        break;
    }
    case 2: {

        break;
    }
    case 3: {
        break;
    }
    case 4: {
        
        break;
    }
    default: {
        std::cout << "Proszę wybrać cyfrę z zakresu wyświetlonego wyżej.";
        break;
    }
    }
        return 0;
}
