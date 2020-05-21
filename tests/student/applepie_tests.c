#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "list.h"
#include "forces.h"
#include "test_util.h"

list_t *make_shape() {
    list_t *shape = list_init(4, free);
    vector_t *v = malloc(sizeof(*v));
    *v = (vector_t){-1, -1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){+1, -1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){+1, +1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){-1, +1};
    list_add(shape, v);
    return shape;
}

// Tests that a mass on a spring oscillates like A cos(sqrt(K / M) * t)
void test_spring_sinusoid() {
    const double M = 10;
    const double K = 2;
    const double A = 3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){A, 0});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);
    for (int i = 0; i < STEPS; i++) {
        assert(vec_isclose(body_get_centroid(mass),
                           (vector_t){A * cos(sqrt(K / M) * i * DT), 0}));
        assert(vec_equal(body_get_centroid(anchor), VEC_ZERO));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

double gravity_potential(double G, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return -G * body_get_mass(body1) * body_get_mass(body2) / sqrt(vec_dot(r, r));
}
double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

// Tests that a gravity force conserves K + U
void test_gravity_energy_conservation() {
    const double M1 = 4.5, M2 = 7.3;
    const double G = 1e3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){10, 20});
    scene_add_body(scene, mass2);
    create_newtonian_gravity(scene, G, mass1, mass2);
    double initial_energy = gravity_potential(G, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        assert(body_get_centroid(mass1).x < body_get_centroid(mass2).x);
        double energy = gravity_potential(G, mass1, mass2) + kinetic_energy(mass1) +
                        kinetic_energy(mass2);
        assert(within(1e-4, energy / initial_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

double spring_potential(double K, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return K * vec_dot(r, r) / 2;
}

// Tests that a spring force conserves K + U
void test_spring_energy_conservation() {
    const double M = 10;
    const double K = 2;
    const double A = 3;
    const double DT = 1e-6;
    const int STEPS = 100000;
    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){A, 0});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);
    double initial_energy = spring_potential(K, mass, anchor);
    for (int i = 0; i < STEPS; i++) {
      double energy = spring_potential(K, mass, anchor) + kinetic_energy(mass);
      assert(within(1e-4, energy / initial_energy, 1));
      scene_tick(scene, DT);
    }
    scene_free(scene);
}

// Tests that a mass with a drag force is slowing down appropriately
void test_drag_location() {
  const double M = 10;
  const double GAMMA = 2;
  const double DT = 1e-6;
  const int STEPS = 100000;
  const int INIT_VEL_X = 5;
  const int INIT_VEL_Y = 5;
  scene_t *scene = scene_init();
  body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
  body_set_velocity(mass, (vector_t) {INIT_VEL_X, INIT_VEL_Y});
  scene_add_body(scene, mass);
  create_drag(scene, GAMMA, mass);
  vector_t prev_vel = body_get_velocity(mass);
  vector_t prev_cent = body_get_centroid(mass);
  assert(vec_isclose(body_get_velocity(mass), (vector_t) {INIT_VEL_X, INIT_VEL_Y}));
  assert(vec_isclose(body_get_centroid(mass), VEC_ZERO));
  scene_tick(scene, DT);
  vector_t new_vel = body_get_velocity(mass);
  vector_t new_cent = body_get_centroid(mass);
  for (int i = 1; i < STEPS; i++) {
    vector_t acc = (vector_t) {-GAMMA * prev_vel.x/ M, -GAMMA * prev_vel.y/ M};
    assert(vec_isclose(new_vel,
                (vector_t){prev_vel.x + acc.x * DT, prev_vel.y + acc.y * DT}));
    assert(vec_isclose(new_cent,
                (vector_t){prev_cent.x + 0.5 * (prev_vel.x + new_vel.x) * DT, prev_cent.y + 0.5 * (prev_vel.y + new_vel.y)* DT}));
    scene_tick(scene, DT);
    prev_vel = new_vel;
    prev_cent = new_cent;
    new_vel = body_get_velocity(mass);
    new_cent = body_get_centroid(mass);
  }
  scene_free(scene);
}

// Tests that a damped harmonic oscillator has correct location
void test_damped_harmonic_oscillator() {
  const double M = 1;
  const double K = 4;
  const double A = 1;
  const double GAMMA = 4;
  const double DT = 1e-3;
  const int STEPS = 5000;
  scene_t *scene = scene_init();
  body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
  body_set_centroid(mass, (vector_t){A, 0});
  scene_add_body(scene, mass);
  body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
  scene_add_body(scene, anchor);
  create_spring(scene, K, mass, anchor);
  create_drag(scene, GAMMA, mass);

  for (int i = 0; i < STEPS; i++) {
    assert(within(1e-2, body_get_centroid(mass).x / (exp(-2*i*DT)*(1 + 2*(i * DT))), 1));

    assert(vec_equal(body_get_centroid(anchor), VEC_ZERO));
    scene_tick(scene, DT);
  }
  scene_free(scene);
}

int main(int argc, char *argv[]) {
    // Run all tests if there are no command-line arguments
    bool all_tests = argc == 1;
    // Read test name from file
    char testname[100];
    if (!all_tests) {
        read_testname(argv[1], testname, sizeof(testname));
    }

    DO_TEST(test_spring_sinusoid)
    DO_TEST(test_gravity_energy_conservation)
    DO_TEST(test_spring_energy_conservation)
    DO_TEST(test_drag_location)
    DO_TEST(test_damped_harmonic_oscillator)

    puts("forces_test PASS");
}
