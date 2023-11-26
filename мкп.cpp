#include <iostream>
#include <cmath>

// Функция для метода простых итераций
double fixedPointIteration(double M, double e, double initial_guess, double tolerance, int max_iterations) {
    double E = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double next_E = M + e * sin(E);
        if (std::abs(next_E - E) < tolerance) {
            return next_E;
        }
        E = next_E;
    }
    return -1; // Возвращает -1 в случае неудачи при нахождении корня
}

// Функция для метода бисекции
double bisectionMethod(double M, double e, double a, double b, double tolerance, int max_iterations) {
    double E_a = a;
    double E_b = b;
    for (int i = 0; i < max_iterations; ++i) {
        double E_mid = (E_a + E_b) / 2;
        double f_mid = E_mid - e * sin(E_mid) - M;
        if (std::abs(f_mid) < tolerance) {
            return E_mid;
        }
        if ((E_a - E_b) / 2 < tolerance) {
            break;
        }
        if ((E_mid - e * sin(E_mid) - M) * (E_a - e * sin(E_a) - M) < 0) {
            E_b = E_mid;
        }
        else {
            E_a = E_mid;
        }
    }
    return -1; // Возвращает -1 в случае неудачи при нахождении корня
}

// Функция для метода золотого сечения
double goldenSectionMethod(double M, double e, double a, double b, double tolerance, int max_iterations) {
    const double golden_ratio = (1 + sqrt(5)) / 2;
    double E_a = a;
    double E_b = b;
    double x1 = E_b - (E_b - E_a) / golden_ratio;
    double x2 = E_a + (E_b - E_a) / golden_ratio;
    for (int i = 0; i < max_iterations; ++i) {
        double f_x1 = x1 - e * sin(x1) - M;
        double f_x2 = x2 - e * sin(x2) - M;
        if (std::abs(E_b - E_a) < tolerance) {
            return (E_a + E_b) / 2;
        }
        if (f_x1 < f_x2) {
            E_b = x2;
            x2 = x1;
            x1 = E_b - (E_b - E_a) / golden_ratio;
        }
        else {
            E_a = x1;
            x1 = x2;
            x2 = E_a + (E_b - E_a) / golden_ratio;
        }
    }
    return -1; // Возвращает -1 в случае неудачи при нахождении корня
}

// Функция для метода Ньютона
double newtonMethod(double M, double e, double initial_guess, double tolerance, int max_iterations) {
    double E = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double f = E - e * sin(E) - M;
        double f_prime = 1 - e * cos(E);
        double delta = f / f_prime;
        E -= delta;
        if (std::abs(delta) < tolerance) {
            return E;
        }
    }
    return -1; // Возвращает -1 в случае неудачи при нахождении корня
}

int main() {
    double M = 0.3; // Средняя аномалия
    double e = 0.8; // Эксцентриситет орбиты

    // Начальное приближение для метода простых итераций
    double initial_guess = M;
    double tolerance = 1e-8; // Точность
    int max_iterations = 1000; // Максимальное количество итераций

    double result_fixed_point = fixedPointIteration(M, e, initial_guess, tolerance, max_iterations);
    double result_bisection = bisectionMethod(M, e, 0, 2 * M, tolerance, max_iterations);
    double result_golden_section = goldenSectionMethod(M, e, 0, 2 * M, tolerance, max_iterations);
    double result_newton = newtonMethod(M, e, initial_guess, tolerance, max_iterations);

    std::cout << "Метод простых итераций: " << result_fixed_point << std::endl;
    std::cout << "Метод бисекции: " << result_bisection << std::endl;
    std::cout << "Метод золотого сечения: " << result_golden_section << std::endl;
    std::cout << "Метод Ньютона: " << result_newton << std::endl;

    return 0;
}