// Translate this file with
//
// g++ -O3 --std=c++11 spacebodies.c -o spacebodies
//
// Run it with
//
// ./spacebodies
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2017 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <cstdio>
#include <limits>

struct Body {
  double mass;
  double velocityX, velocityY, velocityZ;
  double positionX, positionY, positionZ;
  double forceX, forceY, forceZ;
  bool isActive;
};

double t = 0;
double tFinal = 0;
double timeStepSize = 1e-5;

// For dynamic time steps
double smallestDistance = std::numeric_limits<double>::infinity();
int iteration = 2;

int NumberOfBodies = 0;
int NumInactive = 0;

Body *bodies;

std::ofstream videoFile;

// TODO: Remove all debug statements and uncomment Paraview statements

// ----------------------------------------------------
//                    SETUP FUNCTION
// ----------------------------------------------------

void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-2) / 7;

  bodies = new Body[NumberOfBodies];

  int readArgument = 1;

  tFinal = std::stof(argv[readArgument]); readArgument++;

  // Process each body separately
  for (int i=0; i < NumberOfBodies; ++i) {
    // The first three numbers passed into the command line are the x variables
    // representing the position of the particle
    bodies[i].positionX = std::stof(argv[readArgument]); readArgument++;
    bodies[i].positionY = std::stof(argv[readArgument]); readArgument++;
    bodies[i].positionZ = std::stof(argv[readArgument]); readArgument++;

    // The next three numbers passed in represent the velocity vectors
    bodies[i].velocityX = std::stof(argv[readArgument]); readArgument++;
    bodies[i].velocityY = std::stof(argv[readArgument]); readArgument++;
    bodies[i].velocityZ = std::stof(argv[readArgument]); readArgument++;

    bodies[i].mass = std::stof(argv[readArgument]); readArgument++;
    bodies[i].isActive = true;

    if (bodies[i].mass<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
}


// ----------------------------------------------------
//              PARAVIEW RELATED FUNCTIONS
// ----------------------------------------------------

void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}


void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}

/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
void printParaviewSnapshot(int counter) {
  std::stringstream filename;
  filename << "paraview/result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies - NumInactive << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i < NumberOfBodies; ++i) {
    if (!bodies[i].isActive) {
      continue;
    }
    out << bodies[i].positionX
        << " "
        << bodies[i].positionY
        << " "
        << bodies[i].positionZ
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


// ----------------------------------------------------
//                  UPDATE FUNCTIONS
// ----------------------------------------------------

void updatePosition(const double timeStepSize) {
  double accelerationX, accelerationY, accelerationZ;
  double deltaT = (timeStepSize * timeStepSize) / 2;

  for (int i=0; i < NumberOfBodies; ++i) {
    if (bodies[i].isActive) {
      accelerationX = bodies[i].forceX / bodies[i].mass;
      accelerationY = bodies[i].forceY / bodies[i].mass;
      accelerationZ = bodies[i].forceZ / bodies[i].mass;

      bodies[i].positionX = bodies[i].positionX + (timeStepSize * bodies[i].velocityX) + (deltaT * accelerationX);
      bodies[i].positionY = bodies[i].positionY + (timeStepSize * bodies[i].velocityY) + (deltaT * accelerationY);
      bodies[i].positionZ = bodies[i].positionZ + (timeStepSize * bodies[i].velocityZ) + (deltaT * accelerationZ);

      bodies[i].velocityX = bodies[i].velocityX + (timeStepSize * accelerationX);
      bodies[i].velocityY = bodies[i].velocityY + (timeStepSize * accelerationY);
      bodies[i].velocityZ = bodies[i].velocityZ + (timeStepSize * accelerationZ);

      // Set all force to 0 for calculation of next iteration
      bodies[i].forceX = 0;
      bodies[i].forceY = 0;
      bodies[i].forceZ = 0;
    
      // DEBUG
      printf("\nBody %d: %7.8f  %7.8f  %7.8f", i, bodies[i].positionX, bodies[i].positionY, bodies[i].positionZ);
    }
  }
}


void fuseBodies(Body* a, Body* b) {
  double combinedMass = a->mass + b->mass;

  double newVelX = ((a->mass * a->velocityX) + (b->mass * b->velocityX)) / combinedMass;
  double newVelY = ((a->mass * a->velocityY) + (b->mass * b->velocityY)) / combinedMass;
  double newVelZ = ((a->mass * a->velocityZ) + (b->mass * b->velocityZ)) / combinedMass;
  
  // DEBUG
  // printf("\n=====> Old body a: %5.10f, %5.10f, %5.10f, %5.10f", a->mass, a->velocityX, a->velocityY, a->velocityZ);
  // printf("\n=====> Old body b: %5.10f, %5.10f, %5.10f, %5.10f", b->mass, b->velocityX, b->velocityY, b->velocityZ);

  a->mass = combinedMass;
  a->velocityX = newVelX;
  a->velocityY = newVelY;
  a->velocityZ = newVelZ;

  // DEBUG
  printf("\n=====> New combined body : %5.10f, %5.10f, %5.10f, %5.10f", a->mass, a->velocityX, a->velocityY, a->velocityZ);

  b->isActive = false;
  NumInactive += 1;
}


double calculate_distance(Body a, Body b) {
  return sqrt(
    (a.positionX-b.positionX) * (a.positionX-b.positionX) +
    (a.positionY-b.positionY) * (a.positionY-b.positionY) +
    (a.positionZ-b.positionZ) * (a.positionZ-b.positionZ)
  );
}

void addForce(Body* a, Body* b, double distance) {
  double massDistance = a->mass * b->mass / (distance * distance * distance);
  // Calculate the force between body and the others
  a->forceX += (b->positionX-a->positionX) * massDistance ;
  a->forceY += (b->positionY-a->positionY) * massDistance ;
  a->forceZ += (b->positionZ-a->positionZ) * massDistance ;
}

void calcualteEffect(int a_index, int b_index) {
  Body* a = &bodies[a_index];
  Body* b = &bodies[b_index];
  const double distance = calculate_distance(*a, *b);
  // TODO: As distance grows smaller, divide the timestep
  // Keep track of smallest distance
  if (distance < smallestDistance) {
    smallestDistance = distance;
    if (smallestDistance < 0.1) {
      
      // if (smallestDistance < (lastAugmentedDistance / iteration)) {
        double logValue = (1e-6 * log(iteration * distance));
        if (timeStepSize + logValue > 0) {
          timeStepSize += logValue;
          iteration += 1;
        } else if (timeStepSize > 1e-8) {
          // TODO: Figure out how to scale the timestep
          timeStepSize += (1e-6 * log((iteration * 2) * distance));
        }
        // lastAugmentedDistance = smallestDistance / iteration;
        // DEBUG
        std::cout << "\n Time ==== " << timeStepSize << std::endl;
      // }
    }
  }

  // DEBUG
  printf("\nSmallest distance: %5.7f  ------ Distance : %5.7f", smallestDistance, distance); 
  
  // TODO: Change this to 1e-8
  if (distance < 1e-8) {
    // DEBUG
    printf("\n -------------- Bodies %d and %d should collide with distance : %5.7f", a_index, b_index, distance);
    fuseBodies(a, b);
  } else {
    addForce(a, b, distance);
    addForce(b, a, distance);
  }
}

// Part 3: Make the time step change according to how close the bodies are to one another so that the particles don't just pass through each other
void updateBodies() {

  // TODO: Change this timeStepSize;
  // move to global
  const double timeStepSize = 1e-3;
  
  // Step 1.1: All bodies interact and move
  for (int i=0; i < NumberOfBodies; ++i) {
    if (bodies[i].isActive) {
      for (int j=i+1; j < NumberOfBodies; ++j) {
        if (bodies[j].isActive) {
          calcualteEffect(i, j);
        }
      }
    }
  }

  updatePosition(timeStepSize);

  // Increase time
  t += timeStepSize;
}


// ----------------------------------------------------
//                     MAIN METHOD
// ----------------------------------------------------

int main(int argc, char** argv) {
  // Insufficient args
  if (argc==1) {
    std::cerr << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "100.0   0 0 0   1.0 0 0   1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "100.0   0 0 0   1.0 0 0   1.0     0 1.0 0   1.0 0 0   1.0 \t One spiralling around the other one" << std::endl
              << "100.0   3.0 0 0   0 1.0 0   0.4     0 0 0   0 0 0   0.2     2.0 0 0   0 0 0   1.0 \t Three body setup from first lecture" << std::endl
              << std::endl;
    return -1;
  }
  // Mismatched args
  else if ( (argc-2)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  setUp(argc,argv);
  // openParaviewVideoFile();
  // printParaviewSnapshot(0);

  int currentTimeSteps = 0;
  const int plotEveryKthStep = 1;
  while (t<=tFinal) {
    updateBodies();
    currentTimeSteps++;
    if (currentTimeSteps%plotEveryKthStep==0) {
      
      // DEBUG
      // std::cout << "Going into snapshot " << currentTimeSteps/plotEveryKthStep << std::endl;

      // printParaviewSnapshot(currentTimeSteps/plotEveryKthStep);
    }
  }

  // closeParaviewVideoFile();

  return 0;
}