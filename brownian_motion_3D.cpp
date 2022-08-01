#include <random>
#include <Mercury3D.h>
#include <Species/Species.h>
#include <Boundaries/DeletionBoundary.h>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/PeriodicBoundary.h"


using namespace std;

class Brownian3D : public Mercury3D {
public:
    void setupInitialConditions() override {

        const double big_radius = 0.05, small_radius = 0.01;

        // Big Particle
        SphericalParticle p1;
        p1.setSpecies(speciesHandler.getObject(0));
        p1.setRadius(big_radius); // sets particle radius
        // put the big particle in the middle of the box
        p1.setPosition(Vec3D(0.5 * getXMax(), 0.5 * getYMax(), 0.5 * getZMax())); // sets particle position
        p1.setVelocity(Vec3D(0.0, 0.0, 0.0));// sets particle velocity
        particleHandler.copyAndAddObject(p1);

        // Small Particles
        for (int i = 0; i < 2500; i++) {

            // Sample particle positions from uniform distribution but use Gaussian distribution for the velocities
            // This gives us the Maxwell-Boltzmann distribution for the classical gases

            random_device rd;
            mt19937 mt(rd());
            uniform_real_distribution<double> uni_dist(0, 1);
            normal_distribution<double> normal_dist(0, 5);

            double x1 = uni_dist(mt);
            double x2 = uni_dist(mt);
            double x3 = uni_dist(mt);

            // To avoid collision with the big particle, excluded the region below
            if (x1 < getXMax() / 2 + big_radius / 2 && x1 > getXMax() / 2 - big_radius / 2 &&
                x2 < getYMax() / 2 + big_radius / 2 && x2 > getYMax() / 2 - big_radius / 2 &&
                x3 < getZMax() / 2 + big_radius / 2 && x3 > getZMax() / 2 - big_radius / 2) {
                i--;
            } else {
                SphericalParticle p0;
                p0.setSpecies(speciesHandler.getObject(0));
                p0.setRadius(small_radius); // sets particle radius
                p0.setPosition(Vec3D(x1, x2, x3)); // sets particle position
                p0.setVelocity(Vec3D(normal_dist(mt), normal_dist(mt), normal_dist(mt))); // sets particle velocity
                particleHandler.copyAndAddObject(p0);
            }
        }

        // Introduce Periodic Boundary Walls
        PeriodicBoundary b0;
        b0.set(Vec3D(1.0, 0.0, 0.0), 0.0, getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0.0, 1.0, 0.0), 0.0, getYMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0.0, 0.0, 1.0), 0.0, getZMax());
        boundaryHandler.copyAndAddObject(b0);
    }
};


int main(int argc, char *argv[]) {
    Brownian3D problem;

    // Set the general properties of the simulation e.g. name, dimensions, gravity
    problem.setName("brownian_motion_3D");
    problem.setSystemDimensions(3);
    problem.setParticleDimensions(3);
    problem.setGravity(Vec3D(0.0, 0.0, 0.0));
    problem.setXMax(1.0);
    problem.setYMax(1.0);
    problem.setZMax(1.0);
    problem.setTimeMax(1.0);

    // Specify the type of the species and the interaction among them
    LinearViscoelasticFrictionSpecies species;
    species.setDissipation(0.0);
    species.setStiffness(1e6);
    species.setDensity(2000);
    problem.speciesHandler.copyAndAddObject(species);

    // Specify the output file creation frequency
    problem.setSaveCount(10);
    problem.dataFile.setFileType(FileType::ONE_FILE);
    problem.restartFile.setFileType(FileType::ONE_FILE);
    problem.fStatFile.setFileType(FileType::ONE_FILE);
    problem.eneFile.setFileType(FileType::ONE_FILE);
    problem.setWallsWriteVTK(true);
    problem.setParticlesWriteVTK(true);
    logger(INFO, "run number: %", problem.dataFile.getCounter());

    // Set xballs for quick visualization
    problem.setXBallsAdditionalArguments("-solidf -v0 -noborder 4 -cube");

    // Run the DEM solver
    problem.setTimeStep(0.0005 / 50.0);
    problem.solve(argc, argv);

    return 0;
}
