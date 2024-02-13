
#include <random/uniform.h>

typedef float real;

typedef ranlib::MersenneTwister MersenneTwister;
typedef ranlib::independentState independentState;

using ranlib::UniformClosedOpen;

int main(int nargs, char *args[])
{
  UniformClosedOpen<real, MersenneTwister, independentState> r0(0);
  UniformClosedOpen<real, MersenneTwister, independentState> r1(1);
  UniformClosedOpen<real, MersenneTwister, independentState> r2(2);
  UniformClosedOpen<real, MersenneTwister, independentState> r47(47);
  UniformClosedOpen<real, MersenneTwister, independentState> r48(48);
  UniformClosedOpen<real, MersenneTwister, independentState> r49(49);

  r0.seed(1);
  r1.seed(1);
  r2.seed(1);
  r47.seed(1);
  r48.seed(1);
  r49.seed(1);

  for (int i=1; i<10; ++i) {
    std::cout << "i=" << i
              << " r0=" << r0.random()
              << " r0=" << r0.random()
              << " r1=" << r1.random()
              << " r2=" << r2.random()
              << " r47=" << r47.random()
              << " r48=" << r48.random()
              << " r49=" << r49.random()
              << std::endl;
  }

  return 0;
}
