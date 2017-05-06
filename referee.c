#pragma GCC optimize "O3,omit-frame-pointer,inline"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>

int64_t i2 = 0x0000444400004444LL;

#define N_PODS ((int64_t) +4)
#define X_MAX ((int64_t) +16000)
#define Y_MAX ((int64_t) +9000)
#define N_CP_MAX ((int64_t) +8)

#define POD_RADIUS ((int64_t) +400)
#define POD_COLISION_RADIUS ((double) +(POD_RADIUS+POD_RADIUS))
#define POD_COLISION_RADIUS_SQ ((double) +(POD_COLISION_RADIUS*POD_COLISION_RADIUS))

#define CP_RADIUS ((int64_t) +(600))
#define CP_RADIUS_SQ ((double) +CP_RADIUS*CP_RADIUS)
#define MAX_SPEED ((int64_t) +200)

#define FRICTION ((double) +0.85)
#define FRICTION_INT ((int64_t) (+0.85* (1 << 16)))
#define RSPEED ((int64_t) +18)

#define DEPTH ((int64_t) +12)
#define POP ((int64_t) +50)
#define PI ((double) +3.14159265358979323846)

#define SHIELD_CD ((int64_t) +3)

#define ABS(x) ((x)<0 ? -(x) : (x))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// http://files.magusgeek.com/csb/csb_en.html
// https://github.com/sethorizer/csb

struct xy {
    double x;
    double y;
} typedef xy_t;

struct polar {
    double value;
    int64_t angle;
} typedef polar_t;

struct pod {
	xy_t pos;
	xy_t speed;
	xy_t cp;
	xy_t ncp;
	int64_t angle;

	int64_t cpid;
	int64_t ncpid;
	int64_t laps;
	int64_t adv;
	int64_t owner;

	int64_t shielded;
	int64_t shield_cd;
} typedef pod_t;

struct cmd {
	int64_t thrust;
	int64_t angle;
	int64_t shield;
} typedef cmd_t;

struct turn_cmd {
	cmd_t cmds[N_PODS];
} typedef turn_cmd_t;

struct strategy {
	turn_cmd_t cmds[DEPTH];
} typedef strategy_t;

struct population {
	strategy_t indiv[POP];
	double gains[POP];
	int64_t sorted[POP];
} typedef population_t;

struct ctx_t {
	pod_t pods[N_PODS];
	xy_t cps[N_CP_MAX];
	int64_t n_laps;
	int64_t n_cps;
	int64_t turn;
} typedef ctx_t;

ctx_t ctx = {0};
ctx_t ctx_future = {0};
xy_t vect_null = { .x=0.0, .y=0.0 };
xy_t vect_ref = { .x=1.0, .y=0.0 };
xy_t tmp_pos, tmp_speed, tmp1, tmp2, tmp3;

clock_t tstart;
double *current_gain;

strategy_t *strategy;
population_t pop;

double cos_[360] = {1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519, 0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620, 0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146, 0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 0.642788, 0.629320, 0.615661, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 0.500000, 0.484810, 0.469472, 0.453990, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368, 0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069756, 0.052336, 0.034899, 0.017452, 0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104528, -0.121869, -0.139173, -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453990, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798636, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948, -0.891007, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519, -0.951057, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152, -0.809017, -0.798636, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615661, -0.601815, -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472, -0.453990, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568, -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648, -0.156434, -0.139173, -0.121869, -0.104528, -0.087156, -0.069756, -0.052336, -0.034899, -0.017452, -0.000000, 0.017452, 0.034899, 0.052336, 0.069756, 0.087156, 0.104528, 0.121869, 0.139173, 0.156434, 0.173648, 0.190809, 0.207912, 0.224951, 0.241922, 0.258819, 0.275637, 0.292372, 0.309017, 0.325568, 0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453990, 0.469472, 0.484810, 0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615661, 0.629320, 0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913545, 0.920505, 0.927184, 0.933580, 0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848};
double sin_[360] = {0.000000, 0.017452, 0.034899, 0.052336, 0.069756, 0.087156, 0.104528, 0.121869, 0.139173, 0.156434, 0.173648, 0.190809, 0.207912, 0.224951, 0.241922, 0.258819, 0.275637, 0.292372, 0.309017, 0.325568, 0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453990, 0.469472, 0.484810, 0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615661, 0.629320, 0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913545, 0.920505, 0.927184, 0.933580, 0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848, 1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519, 0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620, 0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146, 0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 0.642788, 0.629320, 0.615661, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 0.500000, 0.484810, 0.469472, 0.453990, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368, 0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069756, 0.052336, 0.034899, 0.017452, 0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104528, -0.121869, -0.139173, -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453990, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798636, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948, -0.891007, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519, -0.951057, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152, -0.809017, -0.798636, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615661, -0.601815, -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472, -0.453990, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568, -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648, -0.156434, -0.139173, -0.121869, -0.104528, -0.087156, -0.069756, -0.052336, -0.034899, -0.017452};

void inline add_xy(xy_t *a, xy_t *b) {
	a->x += b->x;
	a->y += b->y;
}

void inline partial_add_xy(xy_t *a, xy_t *b, double f) {
	a->x += b->x*f;
	a->y += b->y*f;
}

void inline sub_xy(xy_t *a, xy_t *b) {
	a->x -= b->x;
	a->y -= b->y;
}

void inline print_xy(xy_t *a) {
    fprintf(stderr, "x:%lf y:%lf ", a->x, a->y);
}

int64_t inline clear_angle(int64_t a) {
	return ((a + 5*180) % 360) - 180;
}

double inline get_dist_sq(xy_t *a, xy_t *b) {
    return (a->x - b->x)*(a->x - b->x) + (a->y - b->y)*(a->y - b->y);
}

double inline get_dist(xy_t *a, xy_t *b) {
    return (double) sqrt(get_dist_sq(a, b));
}

double inline get_angle(xy_t *a, xy_t *b) {
    return (double) (180.0 * atan2(b->y - a->y, b->x - a->x)/PI);
}

double inline get_dotproduct(xy_t *a, xy_t *b) {
    return a->x*b->x + a->y*b->y;
}

double inline get_length(xy_t *a) {
    return get_dist(a, &vect_null);
}

double inline get_angle_vect(xy_t *a, xy_t *b) {
    if (get_length(a) == 0 || get_length(b) == 0) {
        return 0.0;
    }

    double v = get_dotproduct(a, b) / (get_length(a) * get_length(b));
    v = MAX(-1.0, MIN(1.0, v));
    return (double) ((180.0/PI) * acos(v));
}

double inline get_intersect(double dist, int64_t alpha) {
    return (double) (dist*tan(alpha*PI/180.0));
}

polar_t inline get_polar(xy_t *a) {
	polar_t res;
	res.value = get_length(a);
	res.angle = get_angle(&vect_ref, a);
    return res;
}

xy_t inline get_cartesian(polar_t *a) {
	xy_t res;
	res.x = (double) (a->value * cos(PI*a->angle/180));
    res.y = (double) (a->value * sin(PI*a->angle/180));
    return res;
}

void inline do_action(pod_t *pod, cmd_t cmd, int64_t boost) {
	polar_t pol;
	pol.value = 1000.0*MAX_SPEED;
	pol.angle = pod->angle + cmd.angle;

	xy_t target = get_cartesian(&pol);
	add_xy(&target, &(pod->pos));

	printf("%ld %ld ", (int64_t) target.x, (int64_t) target.y);

	if (cmd.shield == 1) {
		printf("SHIELD");
        pod->shield_cd = SHIELD_CD + 1;
	} else if (boost == 1) {
		printf("BOOST");
	} else {
		printf("%ld", cmd.thrust);
	}

	printf(" %ld %ld\n", cmd.angle, cmd.thrust);
}

xy_t inline get_closest(xy_t *a, xy_t *b1, xy_t *b2) {
    // Get clostest point64_t from c that's on the line [a, b]
    xy_t res;
    double da, db, c1, c2, det;

    da = b2->y - b1->y;
    db = b1->x - b2->x;
    c1 = da*b1->x + db*b1->y;
    c2 = -db*a->x + da*a->y;
    det = da*da + db*db;

    if (det != 0) {
        res.x = (da*c1 - db*c2) / det;
        res.y = (da*c2 + db*c1) / det;
    } else {
        // The point is already on the line
        res.x = a->x;
        res.y = a->y;
    }

    return res;
}

double inline get_collision(xy_t *pos1, xy_t *pos2, xy_t *speed1, xy_t *speed2, double radius_sq, double radius) {

    if ((ABS(pos1->x - pos2->x) > 2*radius || ABS(pos1->y - pos2->y) > 2.0*radius) && (ABS(pos1->x + speed1->x - pos2->x - speed2->x) > 2*radius || ABS(pos1->y + speed1->y - pos2->y - speed2->y) > 2.0*radius)) {
    	return -1.0;
    }

    // Square of the distance
    double dist = get_dist_sq(pos1, pos2);

    // We take everything squared to avoid calling sqrt uselessly. It is better for performances

    if (dist < radius_sq) {
        // Objects are already touching each other. We have an immediate collision.
        //fprintf(stderr, "Some shit %lf, %lf\n", dist, radius_sq);
        return -1.0;
    }

    // Optimisation. Objects with the same speed will never collide
    if (speed1->x == speed2->x && speed1->y == speed2->y) {
        return -1.0;
    }

    // We place ourselves in the reference frame of u. u is therefore stationary and is at (0,0)
    double x = pos1->x - pos2->x;
    double y = pos1->y - pos2->y;
    xy_t myp = { .x=x, .y=y };
    double vx = speed1->x - speed2->x;
    double vy = speed1->y - speed2->y;
    xy_t up = {0};
    xy_t toto = { .x=x+vx, .y=y+vy };

    // We look for the closest point to u (which is in (0,0)) on the line described by our speed vector
    xy_t p = get_closest(&up, &myp, &toto);

    // Square of the distance between u and the closest point to u on the line described by our speed vector
    double pdist = get_dist_sq(&up, &p);

    // Square of the distance between us and that point
    double mypdist = get_dist_sq(&myp, &p);

    // If the distance between u and this line is less than the sum of the radii, there might be a collision
    if (pdist < radius_sq) {
     // Our speed on the line
        double length = sqrt(vx*vx + vy*vy);

        // We move along the line to find the point of impact
        double backdist = sqrt(radius_sq - pdist);
        p.x = p.x - backdist * (vx / length);
        p.y = p.y - backdist * (vy / length);

        // If the point is now further away it means we are not going the right way, therefore the collision won't happen
        if (get_dist_sq(&myp, &p) > mypdist) {
            return -1.0;
        }

        pdist = get_dist(&p, &myp);

        // The point of impact is further than what we can travel in one turn
        if (pdist > length) {
            return -1.0;
        }

        if (pdist == 0.0) {
        	return -1.0;
        }

        //fprintf(stderr, "Some other shit %lf, %lf\n", pdist, length);
        // Time needed to reach the impact point
        return pdist / length;
    }

    return -1.0;
}

double inline collision_cp(pod_t* pod) {
	if (pod->cp.x == pod->ncp.x && pod->cp.y == pod->ncp.y) {
		return -1.0;
	}
	return get_collision(&(pod->pos), &(pod->cp), &(pod->speed), &vect_null, CP_RADIUS_SQ, CP_RADIUS);
}

double inline collision_pods(pod_t* pod1, pod_t* pod2) {
	return get_collision(&(pod1->pos), &(pod2->pos), &(pod1->speed), &(pod2->speed), POD_COLISION_RADIUS_SQ, POD_COLISION_RADIUS);
}

void inline bounce_pods(pod_t* pod1, pod_t* pod2) {
    // If a pod has its shield active its mass is 10 otherwise it's 1
    double m1 = pod1->shielded ? 10 : 1;
    double m2 = pod2->shielded ? 10 : 1;
    double mcoeff = (m1 + m2) / (m1 * m2);

    double nx = pod1->pos.x - pod2->pos.x;
    double ny = pod1->pos.y - pod2->pos.y;

    double dvx = pod1->speed.x - pod2->speed.x;
    double dvy = pod1->speed.y - pod2->speed.y;

    // fx and fy are the components of the impact vector. product is just there for optimisation purposes
    double product = nx*dvx + ny*dvy;
    double fx = (nx * product) / (POD_COLISION_RADIUS_SQ * mcoeff);
    double fy = (ny * product) / (POD_COLISION_RADIUS_SQ * mcoeff);

    // We apply the impact vector once
    pod1->speed.x -= fx / m1;
    pod1->speed.y -= fy / m1;
    pod2->speed.x += fx / m2;
    pod2->speed.y += fy / m2;

    // If the norm of the impact vector is less than 120, we normalize it to 120
    double impulse = sqrt(fx*fx + fy*fy);
    if (impulse < 120.0) {
        fx = fx * 120.0 / impulse;
        fy = fy * 120.0 / impulse;
    }

    // We apply the impact vector a second time
    pod1->speed.x -= fx / m1;
    pod1->speed.y -= fy / m1;
    pod2->speed.x += fx / m2;
    pod2->speed.y += fy / m2;
}

void inline init();
void inline reset();
void inline run();

void inline run_full_sim(int64_t target_owner);
void inline run_turn(turn_cmd_t *turn_cmd, int64_t target_owner, int64_t turn);

int64_t main() {
	init();

    while (1) {
    	reset();
    	run();
    }

    return 0;
}

void inline init() {
	int64_t i;
    scanf("%ld", &(ctx.n_laps));
    scanf("%ld", &(ctx.n_cps));

    for (i=0; i < ctx.n_cps; i++) {
        scanf("%lf%lf", &(ctx.cps[i].x), &(ctx.cps[i].y));
    }

    ctx.pods[0].owner = 1;
    ctx.pods[1].owner = 1;
    ctx.pods[2].owner = -1;
    ctx.pods[3].owner = -1;

	for (i=0; i<360; i++) {
		cos_[i] = (double) cos(PI*i/180.0);
		sin_[i] = (double) sin(PI*i/180.0);
    }
}

void inline reset() {
	int64_t i, cpid;
	pod_t *pod;

    for (i=0; i < N_PODS; i++) {
    	//fprintf(stderr, "==== %ld\n", i);
    	pod = &(ctx.pods[i]);
        scanf("%lf%lf%lf%lf%ld%ld", &(pod->pos.x), &(pod->pos.y), &(pod->speed.x), &(pod->speed.y), &(pod->angle), &(cpid));

        if (pod->cpid > 0 && cpid == 0) {
            pod->laps += 1;
        }
        pod->cpid = cpid;
		pod->adv = pod->laps*ctx.n_cps + pod->cpid;
        pod->ncpid = (pod->cpid + 1) % ctx.n_cps;

		//fprintf(stderr, "CP %nd LAP %ld ADV %ld NCP %ld\n", pod->cpid, pod->laps, pod->adv, pod->ncpid);

        pod->cp = ctx.cps[pod->cpid];
        pod->ncp = ctx.cps[pod->ncpid];

        //print_xy(&(pod->cp));
        //print_xy(&(pod->ncp));
        //fprintf(stderr, "\n");
        pod->shield_cd -= 1;
    }
    tstart = clock();
}

void inline run() {
	int64_t i, j;

    if (ctx.turn == 0) {
        printf("%ld %ld BOOST\n", (int64_t) ctx.pods[0].cp.x, (int64_t) ctx.pods[0].cp.y);
        printf("%ld %ld BOOST\n", (int64_t) ctx.pods[0].cp.x, (int64_t) ctx.pods[0].cp.y);
        ctx.turn++;
        return;
    }

    // PUT A GA HERE

	clock_t tend = clock();
    fprintf(stderr, "%f MS\n", (float) (1000*(tend - tstart) / CLOCKS_PER_SEC));

    int64_t winner = 0; // Set the winner as the first individual of the population
    strategy = &(pop.indiv[winner]);
	for (i=0; i < N_PODS; i++) {
		if (ctx.pods[i].owner == -1) {
			continue;
		}
		do_action(&(ctx.pods[i]), strategy->cmds[0].cmds[i], 0);
	}

	ctx.turn += 1;
}

void inline copy_ctx() {
	ctx_future = ctx;
}

void inline run_full_sim(int64_t target_owner) {
    int64_t turn;
    copy_ctx();

    for(turn=0; turn<DEPTH; turn++) {
        run_turn(&(strategy->cmds[turn]), target_owner, turn);
    }
}

void inline run_turn(turn_cmd_t *turn_cmd, int64_t target_owner, int64_t turn) {
	int64_t ipod, jpod, angle, thrust, shield;
	xy_t thrust_vect;
	pod_t *pod;

	// ROTATE AND THRUST
	for (ipod=0; ipod < N_PODS; ipod++) {
		pod = &(ctx_future.pods[ipod]);

		angle = turn_cmd->cmds[ipod].angle;
		thrust = turn_cmd->cmds[ipod].thrust;
        shield = turn_cmd->cmds[ipod].shield;

        angle = (pod->angle + angle + 360) % 360;
        pod->angle = angle;
        pod->shielded = 0;

        if (shield == 1) {
            pod->shielded = 1;
            pod->shield_cd = SHIELD_CD;
            continue;
        }

        if (pod->shield_cd > 0) {
            pod->shield_cd -= 1;
            continue;
        }

        thrust_vect.x = ((double) thrust)*cos_[angle];
        thrust_vect.y = ((double) thrust)*sin_[angle];

        add_xy(&(pod->speed), &thrust_vect);
    }

    // This tracks the time during the turn. The goal is to reach 1.0
    double time = 0.0, fcol_time = 0.0, tmp;
    pod_t *pod1, *pod2;
    int loops = 0;
   	//fprintf(stderr, "AHAH\n");
    // MOVE & COLLISIONS
    while (time < 1.0) {
    	if (loops > 10) {
    		fprintf(stderr, "TIME %lf\n", time);
    	}
        fcol_time = 1.0;
        loops += 1;

        // We look for all the collisions that are going to occur during the turn
        for (ipod=0; ipod < N_PODS; ipod++) {
            // Collision with another pod?
            for (jpod=0; jpod < ipod; jpod++) {
                tmp = collision_pods(&(ctx_future.pods[ipod]), &(ctx_future.pods[jpod]));

            	if (loops > 10) {
            		fprintf(stderr, "DEBUG %lf %lf %ld %ld\n", tmp, fcol_time, ipod, jpod);
            	}

                // If the collision occurs earlier than the one we currently have we keep it
                if (tmp >= 0.0 && tmp + time < 1.0 && tmp < fcol_time) {
                    fcol_time = tmp;
                    pod1 = &(ctx_future.pods[ipod]);
                    pod2 = &(ctx_future.pods[jpod]);
                	if (loops > 10) {
                		fprintf(stderr, "POD COLLISION %lf %lf %ld %ld\n", time, fcol_time, ipod, jpod);
                		fprintf(stderr, "%lf\n", get_dist(&(pod1->pos), &(pod2->pos)));
                	}
                }
            }

            // Collision with another checkpoint?
            // It is unnecessary to check all checkpoints here. We only test the pod's next checkpoint.
            // We could look for the collisions of the pod with all the checkpoints, but if such a collision happens it wouldn't impact the game in any way
            tmp = collision_cp(&(ctx_future.pods[ipod]));

            // If the collision happens earlier than the current one we keep it
            if (tmp >= 0.0 && tmp + time < 1.0 && tmp < fcol_time) {
                fcol_time = tmp;
                pod1 = &(ctx_future.pods[ipod]);
                pod2 = NULL;
                if (loops > 10) {
                	fprintf(stderr, "CP COLLISION %lf %lf %ld\n", time, fcol_time, ipod);
                }
            }
        }

        if (fcol_time >= 1.0) {

            //fprintf(stderr, "END %lf %lf\n", time, fcol_time);

            // No collision, we can move the pods until the end of the turn
            for (ipod=0; ipod < N_PODS; ipod++) {
                pod = &(ctx_future.pods[ipod]);
                if (loops > 10) {
                	fprintf(stderr, "MOVING %ld %lf %lf %lf\n", ipod, pod->speed.x, pod->speed.y, 1.0-time);
                }
                partial_add_xy(&(pod->pos), &(pod->speed), 1.0 - time);
            }

            // End of the turn
            time = 1.0;
        } else {
            // Move the pods to reach the time `t` of the collision
            for (ipod=0; ipod < N_PODS; ipod++) {
                pod = &(ctx_future.pods[ipod]);
                if (loops > 10) {
                	fprintf(stderr, "MOVING %ld %lf %lf %lf\n", ipod, pod->speed.x, pod->speed.y, fcol_time);
                }
                partial_add_xy(&(pod->pos), &(pod->speed), fcol_time);
            }

            if (loops > 10) {
            	fprintf(stderr, "ADVANCING %lf %lf\n", time, fcol_time);
            }

            // Play out the collision
            if (pod2 == NULL) {
		        // New cp reached
	            pod1->cpid = (pod1->cpid + 1) % ctx.n_cps;
	            if (pod1->cpid == 0) {
	                pod1->laps += 1;
	            }

	            pod1->adv = pod1->laps*ctx.n_cps + pod1->cpid;
	            pod1->cp = pod1->ncp;
            } else {
            	bounce_pods(pod1, pod2);
            }

            time += fcol_time;
        }

	    if (loops > 10) {
	    	fprintf(stderr, "END %lf %lf\n", time, fcol_time);
	    }
    }

	// END
	for (ipod=0; ipod < N_PODS; ipod++) {
		pod = &(ctx_future.pods[ipod]);
        pod->speed.x *= FRICTION;
        pod->speed.y *= FRICTION;

        pod->pos.x = round(pod->pos.x);
        pod->pos.y = round(pod->pos.y);
        pod->speed.x = (double) ((int64_t) pod->speed.x);
        pod->speed.y = (double) ((int64_t) pod->speed.y);

        /*if (ipod == 1) {
        	print_xy(&(pod->pos));
        	print_xy(&(pod->speed));
        	fprintf(stderr, "%ld %ld %ld %ld C\n", angle, pod->cpid, pod->laps, pod->adv);
        }*/
	}
    //eval(target_owner, turn);
}
