#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

class particle
{
public:
    particle()
    {
        x.resize(3);
        v.resize(3);
    }
    particle(double xx, double xy, double xz)
    {
        x.resize(3);
        v.resize(3);
        x[0] = xx;
        x[1] = xy;
        x[2] = xz;
    }
    particle(double xx, double xy, double xz, double vx, double vy, double vz)
    {
        x.resize(3);
        v.resize(3);
        x[0] = xx;
        x[1] = xy;
        x[2] = xz;
        v[0] = vx;
        v[1] = vy;
        v[2] = vz;
    }
    ~particle()
    {
    }
    double get_x() { return x[0]; }
    double get_y() { return x[1]; }
    double get_z() { return x[2]; }
    double get_vx() { return v[0]; }
    double get_vy() { return v[1]; }
    double get_vz() { return v[2]; }

private:
    std::vector<double> x;
    std::vector<double> v;
};

class A_class
{
public:
    A_class(double a, double b, double c)
    {
        x = a;
        y = b;
        z = c;
        R = sqrt(x * x + y * y);
        r = sqrt((R - 1) * (R - 1) + z * z);
        e_k.resize(3);
        e_R.resize(3);
        A.resize(3);
        e_k[0] = -y / R;
        e_k[1] = x / R;
        e_k[2] = 0;
        e_R[0] = x / R;
        e_R[1] = y / R;
        e_R[2] = 0;

        A[0] = 0.5 * (0.5 * r * r * e_k[0] / R + 0.5 * z * e_R[0] / R);
        A[1] = 0.5 * (0.5 * r * r * e_k[1] / R + 0.5 * z * e_R[1] / R);
        A[2] = 0.5 * (0.5 * r * r * e_k[2] / R + 0.5 * z * e_R[2] / R - log10(R));
    }
    double get_x() { return A[0]; }
    double get_y() { return A[1]; }
    double get_z() { return A[2]; }

private:
    double x, y, z;
    double R, r;
    std::vector<double> e_k, e_R, A;
};

int main()
{
    // TODO: 模拟结果不正确 nabla A = 0？ nabla phi?
    particle p(1.05, 0, 0, 2.1e-3, 4.3e-4, 0);
    double dt = 0.25; // 步长
    double nt = 1e6;  // 模拟步数
    double m = 1;     // 粒子质量
    double q = 1;     // 电荷
    std::vector<particle> result(nt);
    result[0] = p;
    particle p1(p.get_x() + dt * p.get_vx(), p.get_y() + dt * p.get_vy(), p.get_z() + dt * p.get_vz());
    result[1] = p1;
    for (int i = 1; i < nt - 1; i++)
    {
        double tmp_x, tmp_y, tmp_z;
        auto pl1 = result[i - 1];
        auto pl = result[i];
        A_class Al1(pl1.get_x(), pl1.get_y(), pl1.get_z());
        A_class Al(pl.get_x(), pl.get_y(), pl.get_z());
        tmp_x = q * dt / m * (Al1.get_x() - Al.get_x()) + 2 * pl.get_x() - pl1.get_x();
        tmp_y = q * dt / m * (Al1.get_y() - Al.get_y()) + 2 * pl.get_y() - pl1.get_y();
        tmp_z = q * dt / m * (Al1.get_z() - Al.get_z()) + 2 * pl.get_z() - pl1.get_z();
        result[i + 1] = particle(tmp_x, tmp_y, tmp_z);
    }
    std::ofstream outfile;
    outfile.open("result.txt");
    for (int i = 0; i < nt; i++)
    {

        outfile << result[i].get_x() << "\t" << result[i].get_y() << "\t" << result[i].get_z() << std::endl;
    }
    outfile.close();
    return 0;
}