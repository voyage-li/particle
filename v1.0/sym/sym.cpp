#include <iostream>
#include <vector>
#include <fstream>
int main()
{
    int nt = 1000;
    double dt = 0.05;
    std::vector<double> p(nt), q(nt);
    p[0] = 1;
    q[0] = 0;
    for (int i = 0; i < nt - 1; i++)
    {
        // 谐振子的辛算法模拟
        p[i + 1] = p[i] - dt * q[i];
        q[i + 1] = q[i] + dt * p[i + 1];
    }
    std::ofstream outfile;
    outfile.open("result.txt");
    for (int i = 0; i < nt; i++)
    {
        outfile << q[i] << std::endl;
    }
    outfile.close();
    return 0;
}