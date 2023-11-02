#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#define pi M_PI
int main()
{
    double data[12][3] = {0.1, 1.0152, -4.75613E-15,
                          0.2, 1.06398, -5.51074E-05,
                          0.3, 1.15985, -0.0126204,
                          0.4, 1.28506, -0.066128,
                          0.5, 1.41566, -0.153359,
                          0.6, 1.54571, -0.26411,
                          0.7, 1.67387, -0.392401,
                          0.8, 1.7999, -0.534552,
                          0.9, 1.92387, -0.688109,
                          1.0, 2.0459, -0.85133,
                          1.5, 2.63233, -1.77571,
                          2, 3.18914, -2.8272};
    int id = 7;
    auto k = data[id - 1][0];         // 波数
    auto wr = data[id - 1][1];        // 实部频率
    auto wi = data[id - 1][2];        // 虚部频率
    auto L = 2 * pi / k;              // 模拟空间长度
    auto dt = 0.02;                   // 时间步长
    auto nt = 3000;                   // 模拟的总时间步数
    auto ntout = 200;                 // 输出频率
    auto ng = 32;                     // 离散格点数目
    auto np = 30000;                  // 粒子数量
    auto vb = 1.0;                    // 粒子初始速度
    auto xp1 = 1.0E-2;                // 位置扰动
    auto vp1 = 0.0;                   // 速度扰动
    auto vt = 0.3;                    // 等离子体热速度
    auto wp = 1;                      // 等离子体频率
    auto qm = -1;                     // 荷质比
    auto q = wp * wp / (qm * np / L); // 电荷密度
    auto rho_back = -q * np / L;      // 背景电荷密度
    auto dx = L / ng;                 // 模拟空间离散步长

    std::vector<double> EEk(nt);
    std::vector<double> EEf(nt);
    std::vector<int> t(nt);

    // 初始化位置和速度
    std::vector<double> xp(np);
    std::vector<double> vp(np);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, vt);
    for (int i = 0; i < np; i++)
    {
        xp[i] = (i * L / np);
        vp[i] = distribution(generator) + (1 - 2 * (i % 2)) * vb;
    }

    // 添加扰动
    for (int i = 0; i < np; i++)
    {
        vp[i] += vp1 * cos(k * xp[i]);
        xp[i] += xp1 * cos(k * xp[i]);
    }
    std::vector<int> p(2 * np);
    for (int i = 0; i < np; i++)
    {
        p[i] = i + 1;
        p[i + np] = i + 1;
    }

    for (int it = 0; it < nt; it++)
    {
        for (int i = 0; i < np; i++)
        {
            xp[i] = xp[i] / L + 10.0;
            xp[i] = L * (xp[i] - floor(xp[i]));
        }

        // 绘图
        if (it % (nt / 4) == 0)
        {
            std::cout << it << std::endl;
            std::ofstream outfile;
            outfile.open("./output/" + std::to_string(it / (nt / 4)) + ".txt");

            outfile << it << std::endl;
            outfile << L << std::endl;
            outfile << vt << std::endl;
            outfile << vb << std::endl;
            outfile << dt << std::endl;

            for (int i = 0; i < np; i++)
                outfile << xp[i] << " ";
            outfile << std::endl;
            for (int i = 0; i < np; i++)
                outfile << vp[i] << " ";
            outfile.close();
        }

        // 更新位置
        for (int i = 0; i < np; i++)
        {
            xp[i] = xp[i] + vp[i] * dt;
        }

        std::vector<int> g(np * 2, 0);
        std::vector<int> fraz(np * 2, 0);

        for (int i = 0; i < np; i++)
        {
            g[i] = floor(xp[i] / dx - 0.5) + 1;
            g[i + np] = g[i] + 1;
            fraz[i] = 1 - (double)(fabs(xp[i] / dx - g[i] + 0.5));
            fraz[i + np] = 1 - fraz[i];
        }

        for (int i = 0; i < np * 2; i++)
        {
            if (g[i] < 1)
                g[i] = g[i] + ng;
            else if (g[i] > ng)
                g[i] = g[i] - ng;
        }

        std::vector<std::vector<double>> mat(np, std::vector<double>(ng));
        for (int i = 0; i < 2 * np; i++)
        {
            mat[p[i] - 1][g[i] - 1] = fraz[i];
        }

        std::vector<double> rho(ng, 0.0);
        for (int i = 0; i < ng; i++)
        {
            double sum_mat = 0.0;
            for (int j = 0; j < np; j++)
                sum_mat += mat[j][i];
            rho[i] = q / dx * sum_mat + rho_back;
        }

        std::vector<double> Eg(ng, 0.0);
        double sum_Eg = 0.0;
        for (int j = 0; j < ng - 1; j++)
        {
            Eg[j + 1] = Eg[j] + (rho[j] + rho[j + 1]) * dx / 2;
            sum_Eg += Eg[j + 1];
        }
        Eg[0] = Eg[ng - 1] + rho[ng - 1] * dx;
        sum_Eg += Eg[0];
        for (int j = 0; j < ng; j++)
        {
            Eg[j] = Eg[j] - sum_Eg / ng;
        }

        for (int i = 0; i < np; i++)
        {
            for (int j = 0; j < ng; j++)
            {
                vp[i] += qm * dt * mat[i][j] * Eg[j];
            }
        }

        double eek_it = 0.0;
        double eef_it = 0.0;
        for (int i = 0; i < np; i++)
            eek_it += 0.5 * (double)fabs(q) * vp[i] * vp[i];
        for (int i = 0; i < ng; i++)
            eef_it += 0.5 * dx * Eg[i] * Eg[i];

        EEk[it] = eek_it;
        EEf[it] = eef_it;
        t[it] = it * dt;
    }
    return 0;
}