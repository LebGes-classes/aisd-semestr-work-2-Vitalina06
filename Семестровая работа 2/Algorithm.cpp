#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <cmath>
#include <limits>

using namespace std;

struct Point {
    double x, y;
    
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    // Перегрузка оператора для сортировки точек
    bool operator < (const Point& p) const {
        return x < p.x || (x == p.x && y < p.y);
    }
    
    // Перегрузка оператора для сравнения точек
    bool operator == (const Point& p) const {
        return x == p.x && y == p.y;
    }
};

// Функция для вычисления векторного произведения (для определения ориентации)
double cross_product(const Point& O, const Point& A, const Point& B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Алгоритм Джарвиса (используется внутри алгоритма Чана)
vector<Point> jarvis_march(const vector<Point>& points) {
    int n = points.size();
    if (n <= 3) return points;
    
    vector<Point> hull;
    
    // Находим самую левую точку
    int leftmost = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].x < points[leftmost].x || 
           (points[i].x == points[leftmost].x && points[i].y < points[leftmost].y)) {
            leftmost = i;
        }
    }
    
    int p = leftmost, q;
    do {
        hull.push_back(points[p]);
        
        q = (p + 1) % n;
        for (int i = 0; i < n; i++) {
            // Если i находится "левее" pq, то обновляем q
            if (cross_product(points[p], points[i], points[q]) < 0) {
                q = i;
            }
        }
        
        p = q;
    } while (p != leftmost);
    
    return hull;
}

// Алгоритм Грэхема (используется внутри алгоритма Чана)
vector<Point> graham_scan(vector<Point>& points) {
    int n = points.size();
    if (n <= 3) return points;
    
    // Находим самую нижнюю точку (и самую левую, если таких несколько)
    int min_idx = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].y < points[min_idx].y || 
           (points[i].y == points[min_idx].y && points[i].x < points[min_idx].x)) {
            min_idx = i;
        }
    }
    
    // Помещаем самую нижнюю точку в начало
    swap(points[0], points[min_idx]);
    
    // Сортируем точки по полярному углу относительно points[0]
    Point p0 = points[0];
    sort(points.begin() + 1, points.end(), [p0](const Point& a, const Point& b) {
        double cross = cross_product(p0, a, b);
        if (cross == 0) {
            // Если точки коллинеарны, сортируем по расстоянию
            return (a.x - p0.x)*(a.x - p0.x) + (a.y - p0.y)*(a.y - p0.y) < 
                   (b.x - p0.x)*(b.x - p0.x) + (b.y - p0.y)*(b.y - p0.y);
        }
        return cross > 0;
    });
    
    // Удаляем коллинеарные точки (оставляем самую дальнюю)
    vector<Point> unique_points;
    unique_points.push_back(points[0]);
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && cross_product(p0, points[i], points[i+1]) == 0) {
            i++;
        }
        unique_points.push_back(points[i]);
    }
    
    if (unique_points.size() < 3) return unique_points;
    
    // Строим выпуклую оболочку
    stack<Point> hull_stack;
    hull_stack.push(unique_points[0]);
    hull_stack.push(unique_points[1]);
    hull_stack.push(unique_points[2]);
    
    for (int i = 3; i < unique_points.size(); i++) {
        while (hull_stack.size() >= 2) {
            Point top = hull_stack.top();
            hull_stack.pop();
            Point next_top = hull_stack.top();
            
            if (cross_product(next_top, top, unique_points[i]) > 0) {
                hull_stack.push(top);
                break;
            }
        }
        hull_stack.push(unique_points[i]);
    }
    
    // Переносим результаты из стека в вектор
    vector<Point> hull;
    while (!hull_stack.empty()) {
        hull.push_back(hull_stack.top());
        hull_stack.pop();
    }
    
    // Разворачиваем, так как стек возвращает точки в обратном порядке
    reverse(hull.begin(), hull.end());
    
    return hull;
}

// Алгоритм Чана для построения выпуклой оболочки
vector<Point> chans_algorithm(const vector<Point>& points) {
    int n = points.size();
    if (n <= 5) {
        // Для небольшого количества точек используем алгоритм Джарвиса
        return jarvis_march(points);
    }
    
    // Выбираем параметр m (количество подмножеств)
    int m = min(n, (int)ceil(exp(log(n))));
    
    // Разбиваем точки на m подмножеств по n/m точек в каждом
    vector<vector<Point>> subsets(m);
    for (int i = 0; i < n; i++) {
        subsets[i % m].push_back(points[i]);
    }
    
    // Для каждого подмножества строим выпуклую оболочку алгоритмом Грэхема
    vector<vector<Point>> convex_hulls(m);
    for (int i = 0; i < m; i++) {
        convex_hulls[i] = graham_scan(subsets[i]);
    }
    
    // Объединяем все точки из всех выпуклых оболочек
    vector<Point> all_hull_points;
    for (const auto& hull : convex_hulls) {
        for (const auto& p : hull) {
            all_hull_points.push_back(p);
        }
    }
    
    // Применяем алгоритм Джарвиса к объединенным точкам
    return jarvis_march(all_hull_points);
}

int main() {
    // Пример использования
    vector<Point> points = {
        {0, 0}, {1, 1}, {2, 2}, {3, 1}, {4, 0},
        {3, -1}, {2, -2}, {1, -1}, {0, -2}, {1, 0}
    };
    
    vector<Point> convex_hull = chans_algorithm(points);
    
    cout << "Выпуклая оболочка содержит следующие точки:\n";
    for (const auto& p : convex_hull) {
        cout << "(" << p.x << ", " << p.y << ")\n";
    }
    
    return 0;
}