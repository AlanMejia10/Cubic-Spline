#include <iostream>
#include <cmath>
#include <vector>
#include <deque>
#include <algorithm>
#include <armadillo>

// Prototipos de funciones

bool OrdernarPorValor(const std::pair<float, float> &a, const std::pair<float, float> &b) {
    return (a.first < b.first);
}

void SplineCubico();
void clearScreen();
void SetupTable();
void SortData();
std::vector<arma::rowvec> Spline(std::vector<float>& anchuraH, std::vector<float>& fx);
void LeerTabla();
void ImprimirTabla();
void ModificarPosicion(int position);

// Variable global
std::vector<std::pair<float, float>> data;

int main(){
    std::cout << "Spline cubico" << std::endl;
    // -------- Aqui pon lo que debas poner ---------
    SplineCubico();
}

void SplineCubico() {
    std::vector<float> anchuraH;
    std::vector<float> fx;

    clearScreen();
    SetupTable();
    SortData();

    //Matrix matCoeficientes = Spline(anchuraH, fx);
    std::vector<arma::rowvec> matCoeficientes = Spline(anchuraH, fx);

    for (int i = 0; i < (int)data.size() - 1; i++) {
        std::cout << "g" << i << "(x)=" << std::showpos << matCoeficientes[i](0) << "x^3 " << matCoeficientes[i](1) << "x^2 " << matCoeficientes[i](2) << "1x " << matCoeficientes[i](3) << std::endl;
        std::cout << data[i].first << " <= x <= " << data[i + 1].first << "\n"<< std::endl;
        std::cout << std::noshowpos;
    }
}


void clearScreen() {
#ifdef _WIN32
    system("cls");
#elif defined(unix) || defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
    system("clear");
#else
#error "OS not supported."
#endif
}

void SetupTable() {
    char opc;

    LeerTabla();
    clearScreen();
    ImprimirTabla();
    std::cout << "Los datos de la tabla son correctos: (s/n): ";
    std::cin >> opc;

    if (opc == 's')
        return;

    while (opc != 's') {
        int position;
        std::cout << "Indica la posicion que se modificara: ";
        std::cin >> position;
        ModificarPosicion(position);
        ImprimirTabla();
        std::cout << "Los datos de la tabla son correctos: (s/n): ";
        std::cin >> opc;
    }
}

void SortData() {
    sort(data.begin(), data.end(), OrdernarPorValor);
}

std::vector<arma::rowvec> Spline(std::vector<float>& anchuraH, std::vector<float>& fx) {
    for (int i = 0; i < (int)data.size() - 1; i++) {
        float deltaH = data[i + 1].first - data[i].first;
        float func = (data[i + 1].second - data[i].second) / deltaH;
        anchuraH.push_back(deltaH);
        fx.push_back(func);
    }

    arma::mat matA(data.size() - 2, data.size() - 2);

    for (int i = 0; i < (int)matA.n_rows; ++i) {
        for (int j = 0; j < (int)matA.n_cols; ++j) {
            if (j == i - 1)
                matA(i,j) = anchuraH[i];
            else if (i == j)
                matA(i,j) = 2 * (anchuraH[i] + anchuraH[i + 1]);
            else if (j == i + 1)
                matA(i,j) = anchuraH[i + 1];
            else
                matA(i,j) = 0;
        }
    }

    arma::mat matB(data.size() - 2, 1);

    for (int i = 0; i < (int)matB.n_rows; i++)
        matB(i,0) = 6 * (fx[i + 1] - fx[i]);

    arma::mat matInv = arma::inv(matA);
    arma::mat matResult = matInv * matB;

    std::deque<float> sCoeficientes;
    for (int i = 0; i < (int)matResult.n_rows; i++) {
        sCoeficientes.push_back(matResult(i,0));
    }
    sCoeficientes.push_front(0);
    sCoeficientes.push_back(0);


    // ----- Calculando valores de ai, bi, ci y di
    arma::mat matCoeficientes(anchuraH.size(), 4);
    for (int i = 0; i < (int)matCoeficientes.n_rows; i++) {
        for (int j = 0; j < (int)matCoeficientes.n_cols; j++) {
            if (j == 0)
                matCoeficientes(i,j) = (sCoeficientes[i + 1] - sCoeficientes[i]) / (6 * anchuraH[i]);
            else if (j == 1)
                matCoeficientes(i,j) = sCoeficientes[i] / 2.0;
            else if (j == 2)
                matCoeficientes(i,j) = ((data[i + 1].second - data[i].second) / anchuraH[i]) - anchuraH[i] * ((sCoeficientes[i + 1] + 2 * sCoeficientes[i]) / 6.0);
            else if (j == 3)
                matCoeficientes(i,j) = data[i].second;
        }
    }

    // ---------- calcular producto notable

    std::vector<arma::rowvec> coeficientesSpline;
    for(int i = 0; i < (int) matCoeficientes.n_rows; i++) {
        arma::rowvec coef(4);
        float valorB = data[i].first * -1;
        coef(0) = matCoeficientes(i, 0);
        coef(1) = matCoeficientes(i, 0) * 3 * valorB + matCoeficientes(i, 1);
        coef(2) = matCoeficientes(i, 0) * 3 * pow(valorB, 2) + matCoeficientes(i, 1) * 2 * valorB + matCoeficientes(i, 2);
        coef(3) = matCoeficientes(i, 0) * pow(valorB, 3) + matCoeficientes(i, 1) * pow(valorB, 2) + matCoeficientes(i, 2) * valorB + matCoeficientes(i, 3);
        coeficientesSpline.push_back(coef);
    }

    return coeficientesSpline;
}

void LeerTabla() {
    int numPuntos;
    float x, fx;
    std::cout << "Ingresa el numero de puntos: ";
    std::cin >> numPuntos;

    for (int i = 0; i < numPuntos; ++i) {
        std::cout << "Ingresa el valor de x" << i << ": ";
        std::cin >> x;
        std::cout << "Ingresa el valor de f" << i << ": ";
        std::cin >> fx;
        data.push_back(std::make_pair(x, fx));
    }
}

void ImprimirTabla() {
    std::cout << "i \tx \tf(x) " << std::endl;
    for (int i = 0; i < (int)data.size(); i++)
        std::cout << i << "\t" << data[i].first << "    " << data[i].second << std::endl;
}

void ModificarPosicion(int position) {
    float x, fx;
    std::cout << "Ingresa el valor de x" << position << ": ";
    std::cin >> x;
    std::cout << "Ingresa el valor de f" << position << ": ";
    std::cin >> fx;
    data.at(position) = std::make_pair(x, fx);
}
