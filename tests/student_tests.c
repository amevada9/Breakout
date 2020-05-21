// // TODO: IMPLEMENT YOUR TESTS IN THIS FILE
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "forces.h"
#include "test_util.h"
#include "vector.h"

const vector_t WINDOW = {1000, 500};


  // list_t *make_body() {
  //     // not really sure what the veolcities and all should be
  //     vector_t center = {WINDOW.x / 2, WINDOW.y /2};
  //     list_t *shape = list_init(4, free);
  //     vector_t *v = malloc(sizeof(*v));
  //     *v = (vector_t){center.x-1, center.y-1};
  //     list_add(shape, v);
  //     v = malloc(sizeof(*v));
  //     *v = (vector_t){center.x+1, center.y-1};
  //     list_add(shape, v);
  //     v = malloc(sizeof(*v));
  //     *v = (vector_t){center.x+1, center.y+1};
  //     list_add(shape, v);
  //     v = malloc(sizeof(*v));
  //     *v = (vector_t){center.x-1, center.y+1};
  //     list_add(shape, v);
  //     return shape;
  //   }


    list_t *make_body() {
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



  //simple kinematics with a simple force
  void test_orbit(){
    const double G_CONST = 6.67e-11;
    const double MASS_1 = 1;
    const double MASS_2 = 100000;
    const int STEPS = 1000000;
    //const double DT = 1e-6;


    scene_t *scene = scene_init();
    body_t *orbit = body_init(make_body(),  MASS_1, (rgb_color_t){0, 0, 0});
    body_set_centroid(orbit, (vector_t) {0,10}); // might be a prole,?
    body_t *mass = body_init(make_body(),  MASS_2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t) {0, 0});
    vector_t orbit_vel = {100, 0};
    body_set_velocity(orbit, orbit_vel);
    //vector_t position2 = {1000 ,0};
    scene_add_body(scene, orbit);
    scene_add_body(scene, mass);
    create_newtonian_gravity(scene, G_CONST, orbit, mass);
    double init_rad = 10.0;
      for (int i = 0; i < STEPS; i++) {
        vector_t centroid_body_1 = body_get_centroid(orbit);
        vector_t centroid_body_2 = body_get_centroid(mass);
        vector_t difference_centroids = vec_subtract(centroid_body_2, centroid_body_1);
        double distance_bodies = vec_magnitude(difference_centroids);

        assert(within(1e-4, distance_bodies / init_rad, 1));

        // size_t mass_body_1 = body_get_mass(body1);
        // size_t mass_body_2 = body_get_mass(body2);
        // vector_t velo1 = body_get_velocity(body1);
        // vector_t velo2 = body_get_velocity(body2);
        // double forceMutual = (G_CONST * mass_body_1 * mass_body_2 )/ (distance_bodies * distance_bodies);
        // double acc1 = forceMutual / mass_body_1;
        // double acc2 = forceMutual / mass_body_2;
        // vector_t vec_change_1 = {velo1.x * DT + 0.5 + acc1 * DT * DT, 0};
        // position1 = vec_add(position1, vec_change_1);
        // vector_t vec_change_2 = {velo2.x * DT + 0.5 + acc2 * DT * DT, 0};
        // position2 = vec_add(position2, vec_change_2);
        // scene_tick(scene, DT);
        // assert(vec_isclose(body_get_centroid(body1), position1));
        // assert(vec_isclose(body_get_centroid(body2), position2));
      }
    scene_free(scene);
  }

  //energy concseravtion in spring (with drag if wanted)
  void test_spring_energy_conservation(){
    const double MASS_1 = 10;
    const double MASS_2= 20;
    const double K = 2.0;
    const double A = 3.0;
    const double DT = 1e-6;
    const int STEPS = 10000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_body(), MASS_1, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass1, (vector_t){A, 0});
    scene_add_body(scene, mass1);
    //infinity instead of big number
    body_t *mass2 = body_init(make_body(), MASS_2, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass2);
    create_spring(scene, K, mass1, mass2);

    double total_energy = .5 * K * A * A; // think I'm interpeting this right
    for (int i = 0; i < STEPS; i++) {
      vector_t centroid1 = body_get_centroid(mass1);
      vector_t centroid2 = body_get_centroid(mass2);
      vector_t difference = vec_subtract(centroid1, centroid2);
      double distance = vec_magnitude (difference);
      double energy = .5 * K * distance * distance;
      // printf(" calc energy: %f\n",energy  );
      // printf(" total energy: %f\n",total_energy  );

        assert(within(1e-4, energy / total_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}


  // drag with simple movement (unitl it stops)
  void test_drag(){
    const double DRAG = .75;
    const double MASS = 2;
    const int STEPS = 1000000;
    const vector_t INIT_VEL = {5000, 0};
    const double DT = 1e-6;

    scene_t *scene = scene_init();
    //vector_t center = vec_init(WINDOW.x/2, WINDOW.y/2);
    body_t *body = body_init(make_body(), MASS, (rgb_color_t){0, 0, 0});
    body_set_centroid(body, (vector_t){WINDOW.x/2, WINDOW.y/2});
    body_set_velocity(body, (vector_t){INIT_VEL.x, INIT_VEL.y}); // not really sure if this vec_init is necessary
    scene_add_body(scene, body);
    create_drag(scene, DRAG, body);
    double x = WINDOW.x /2;
    vector_t velocity = INIT_VEL;
    for (int i = 0; i < STEPS; i++) {
        scene_tick(scene, DT);
      // calculating new drag force
        double force = - DRAG * velocity.x;
        double drag_acc = force / body_get_mass(body);


        vector_t updated_v = {velocity.x + drag_acc * DT, 0};

        // updating expected centroid position
        double updated_x = x + velocity.x * DT + .5 * drag_acc * DT * DT;
        x = updated_x;
        velocity = updated_v;

        // updating the resultant velocity
        vector_t new_pos = {x, WINDOW.y/2};
        assert(vec_isclose(body_get_centroid(body), new_pos));

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

      // DO_TEST(test_simple_gravity)
      // DO_TEST(test_spring_energy_conservation)

      DO_TEST(test_orbit)
      DO_TEST(test_drag)
      DO_TEST(test_spring_energy_conservation)

      puts("forces_test PASS");
  }
