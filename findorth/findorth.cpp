#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const double pi = 4.0 * atan(1.0);

static double rel_err(double approx, double real)
{
    return fabs(approx - real) / real;
}

static int gcd(int a, int b)
{
    if(b == 0)
        return a;
    else
        return gcd(b, a % b);
}

static double log2(double x)
{
    return log(x) / log(2.0);
}

struct Cost
{
    int adds;
    int shifts;

    Cost() : adds(0), shifts(0) {}
    Cost(int a, int s) : adds(a), shifts(s) {}

    int sum_cost() const { return adds + shifts; }
};

Cost operator *(int s, const Cost &a)
{
    return Cost(s * a.adds, s * a.shifts);
}

Cost operator +(const Cost &a, const Cost &b)
{
    return Cost(a.adds + b.adds, a.shifts + b.shifts);
}

static Cost get_mult_cost(int val, int scale_log2)
{
    Cost c;
    int scale = 1 << scale_log2;
    int hadbits = 0;

    // if the "scale" bit is set but the one below it isn't, it's cheaper
    // to force it to be treated as an isolated bit, rather than treating
    // the whole thing is an interrupted run, because the scale bit is the
    // one bit we can produce without a shift.
    // 
    // in particular, an input like (binary):
    //   1.011
    // will be considered as an interrupted 1-run:
    //   (10.000 - 0.001) - 0.100
    // instead of the cheaper (1 shift less)
    //   1.000 + (0.100 - 0.001)
    //
    // hence, special-case this here.
    if ((val & scale) != 0 && (val & (scale >> 1)) == 0) {
        hadbits = 1;
        val &= ~scale;
    }

    // this is the bits we have to add/subtract stuff at
    int events = val ^ (val << 1);

    while (events) {
        int this_bit = events & -events;
        events &= events - 1;

        // is this a single-bit interuption in a run?
        if (events & (this_bit << 1)) {
            // yes it is. cheaper to treat this as a single bit,
            // not a 1-bit run; this boils down to "mask out the
            // trailing edge".
            events &= events - 1;
        }

        c.adds += hadbits;
        c.shifts += this_bit != scale;
        hadbits = 1;
    }

    return c;
}

// Find the cheapest way to compute the rotation. See below.
static const char *best_rot_variant(int c, int s, int scale_log2, Cost *out_cost)
{
    Cost cost_c = get_mult_cost(c, scale_log2);
    Cost cost_s = get_mult_cost(s, scale_log2);
    Cost cost_sum = get_mult_cost(c + s, scale_log2);
    Cost cost_diff = get_mult_cost(abs(c - s), scale_log2);

    // direct matrix multiply variant
    const char *var = "direct";
    *out_cost = 2*cost_c + 2*cost_s + Cost(2, 0); // 4 muls plus 2 adds.

    // c*(x+y) variant
    Cost cost_cxpy = cost_c + cost_sum + cost_diff + Cost(3, 0); // 3 muls plus 3 adds.
    if (cost_cxpy.sum_cost() < out_cost->sum_cost()) {
        var = "c*(x+y)";
        *out_cost = cost_cxpy;
    }

    // s*(x-y) variant
    Cost cost_sxmy = cost_s + cost_diff + cost_sum + Cost(3, 0); // 3 muls plus 3 adds
    if (cost_sxmy.sum_cost() < out_cost->sum_cost()) {
        var = "s*(x-y)";
        *out_cost = cost_sxmy;
    }

    return var;
}

// This function finds a pair of approximations
// [c1 -s1]       [c2 -s2]
// [s1  c1]  and  [s2  c2]
// to given planar rotations subject to the following constraints:
// 1. c1, s1, c2, s2 are all integer
// 2. c1*c1 + s1*s1 == c2*c2 + s2*s2, i.e. the approximations have the same norm.
// This is done while simultaneously trying to minimize:
// 1. The relative errors of the approximations
// 2. The size of the integers in question (smaller is better)
// 3. The cost of computing this approximation when no integer multiplies are available.
//    Computing the minimum cost of an integer multiply is fairly tricky in and of itself.
//    Here, an added complication is that this is still a (scaled) rotation matrix, and
//    the rotation
//      x' = c*x - s*y
//      y' = s*x + c*y
//    can be equivalently factored as either
//      t = c*(x+y)
//      x' = t - (c+s)*y
//      y' = t + (s-c)*x
//    or
//      t = s*(x-y)
//      x' = t + (c-s)*x
//      y' = t + (c+s)*y
//    which has less integer multiplies and might be cheaper to compute, depending on how the
//    factors come out. Oh, and there might be a way to share computations between the products,
//    but this function doesn't consider that.
static void find_rot_pair(int r0, int r1)
{
    double c1 = cos((double)r0 * pi / 16.0);
    double s1 = sin((double)r0 * pi / 16.0);
    double c3 = cos((double)r1 * pi / 16.0);
    double s3 = sin((double)r1 * pi / 16.0);

    double s1_over_c1 = s1 / c1;
    double c3_over_c1 = c3 / c1;
    double s3_over_c1 = s3 / c1;

    for (int ic1=1; ic1<256; ic1++) {
        int vals[3];
        vals[0] = (int)(ic1 * s1_over_c1);
        vals[1] = (int)(ic1 * c3_over_c1);
        vals[2] = (int)(ic1 * s3_over_c1);

        double s = c1 / ic1;

        for (int round=0; round<8; round++) {
            int is1 = vals[0] + ((round >> 0) & 1);
            int ic3 = vals[1] + ((round >> 1) & 1);
            int is3 = vals[2] + ((round >> 2) & 1);

            // only consider minimal fractions
            if (gcd(ic1, gcd(is1, gcd(ic3, is3))) != 1)
                continue;

            if (ic1*ic1 + is1*is1 == ic3*ic3 + is3*is3) {
                Cost rot1c, rot1cb;
                Cost rot3c, rot3cb;
                int scale_log2 = (int) log2(ic1);

                const char *var1 = best_rot_variant(ic1, is1, scale_log2, &rot1c);
                const char *var3 = best_rot_variant(ic3, is3, scale_log2, &rot3c);
                Cost sum_cost = rot1c + rot3c;

                const char *var1b = best_rot_variant(ic1, is1, scale_log2 + 1, &rot1cb);
                const char *var3b = best_rot_variant(ic3, is3, scale_log2 + 1, &rot3cb);
                Cost sum_costb = rot1cb + rot3cb;
                if (sum_costb.sum_cost() < sum_cost.sum_cost()) {
                    scale_log2++;
                    var1 = var1b;
                    var3 = var3b;
                    rot1c = rot1cb;
                    rot3c = rot3cb;
                    sum_cost = sum_costb;
                }

                printf("candidate: (%d,%d), (%d,%d) norm: %.2f  rel err: %.4f %.4f %.4f\n", ic1, is1, ic3, is3,
                        hypot(ic1, is1), rel_err(is1 * s, s1), rel_err(ic3 * s, c3), rel_err(is3 * s, s3));

                printf("  sum cost: %dA %dS  scale: %d\n", sum_cost.adds, sum_cost.shifts, 1 << scale_log2);
                printf("  rot 1: %s %dA %dS\n  rot 3: %s %dA %dS\n",
                        var1, rot1c.adds, rot1c.shifts,
                        var3, rot3c.adds, rot3c.shifts);
            }
        }
    }
}

static void find_rot(int r)
{
    double c2 = cos((double)r * pi / 16.0);
    double s2 = sin((double)r * pi / 16.0);

    double s2_over_c2 = s2 / c2;

    for (int ic2=1; ic2<128; ic2++) {
        int val = (int)(ic2 * s2_over_c2);

        double s = c2/ic2;

        for (int round=0; round<2; round++) {
            int is2 = val + round;
            if (gcd(ic2, is2) != 1) // only consider solutions in minimal terms
                continue;

            Cost rot2c, rot2cb;
            int scale_log2 = (int) log2(ic2);

            const char *var2 = best_rot_variant(ic2, is2, scale_log2, &rot2c);
            const char *var2b = best_rot_variant(ic2, is2, scale_log2 + 1, &rot2cb);
            
            if (rot2cb.sum_cost() < rot2c.sum_cost()) {
                scale_log2++;
                var2 = var2b;
                rot2c = rot2cb;
            }

            double err = rel_err(is2 * s, s2);
            if (err < 0.06 && rot2c.sum_cost() <= 11) {
                printf("candidate: (%d,%d) norm: %.2f  rel err: %.4f\n", ic2, is2, hypot(ic2, is2), err);
                printf("  rot: %s %dA %dS scale=%d\n", var2, rot2c.adds, rot2c.shifts, 1 << scale_log2);
            }
        }
    }
}

int main()
{
    printf("rot2:\n");
    find_rot(2);

    printf("\nrot1/3 pair:\n");
    find_rot_pair(1, 3);

    return 0;
}
