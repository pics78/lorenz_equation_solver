#include <fstream>
#include <iomanip>
#define DATA_FILE "lorenz_eq.dat"

// Lorenz equation
double dx_dt(double p, double x, double y) { return -p*(x - y); }
double dy_dt(double r, double x, double y, double z) { return -x*z + r*x - y; }
double dz_dt(double b, double x, double y, double z) { return x*y - b*z; }

int main() {

  // 3 constants
  const double p = 10.0;
  const double r = 28.0;
  const double b = 8/3;

  // initial position
  const double x0 = -1.0;
  const double y0 = 2.0;
  const double z0 = 1.0;

  // calculation condition
  const double dt    = 0.001;
  const double t_max = 20.0;

  const int timestep = t_max / dt;

  double x = x0;
  double y = y0;
  double z = z0;
  double t = 0.0;

  std::ofstream outputfile;
  outputfile.open(DATA_FILE);
  outputfile << std::fixed << std::setprecision(7)
             << "### lorenz equation - numerical calculation result : (x, y, z) ###\n"
             << "# const : (p, r, b) = (" << p << ", " << r << ", " << b << ")\n\n"
             << "# t = " << t << "\n"
             << x0 << "\t" << y0 << "\t" << z0 << "\n";

  double frac_6 = 1.0 / 6.0;

  for (int i=0; i<timestep; i++) {

    t += dt;

    // RK4 method
    double k1 = dt * dx_dt(p, x, y);          double k1x = x + 0.5*k1;
    double l1 = dt * dy_dt(r, x, y, z);       double l1y = y + 0.5*l1;
    double m1 = dt * dz_dt(b, x, y, z);       double m1z = z + 0.5*m1;

    double k2 = dt * dx_dt(p, k1x, l1y);      double k2x = x + 0.5*k2;
    double l2 = dt * dy_dt(r, k1x, l1y, m1z); double l2y = y + 0.5*l2;
    double m2 = dt * dz_dt(b, k1x, l1y, m1z); double m2z = z + 0.5*m2;

    double k3 = dt * dx_dt(p, k2x, l2y);      double k3x = x + k3;
    double l3 = dt * dy_dt(r, k2x, l2y, m2z); double l3y = y + l3;
    double m3 = dt * dz_dt(b, k2x, l2y, m2z); double m3z = z + m3;

    double k4 = dt * dx_dt(p, k3x, l3y);
    double l4 = dt * dy_dt(r, k3x, l3y, m3z);
    double m4 = dt * dz_dt(b, k3x, l3y, m3z);

    // update (x, y, z)
    x = x + (k1 + 2.0*k2 + 2.0*k3 + k4) * frac_6;
    y = y + (l1 + 2.0*l2 + 2.0*l3 + l4) * frac_6;
    z = z + (m1 + 2.0*m2 + 2.0*m3 + m4) * frac_6;

    outputfile << "# t = " << t << "\n"
               << x << "\t" << y << "\t" << z << "\n";
  }

  outputfile.close();
  return 0;
}
