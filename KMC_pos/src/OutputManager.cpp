#include "OutputManager.h"
#include <fstream>

void OutputManager::registerTOF(const double time) {
  std::ofstream outFile;
  outFile.open(tofFile, std::fstream::app);

  if (outFile.is_open()) {
    outFile << boost::format("%10.10f ") % time;
  } else {
    std::cout << "Could not open output file" << std::endl;
  }
}

void OutputManager::endTOF() {
  std::ofstream outFile;
  outFile.open(tofFile, std::fstream::app);
  if (outFile.is_open()) {
    outFile << std::endl;
  } else {
    std::cout << "Could not open output file" << std::endl;
  }
}

void OutputManager::registerParticlePositions(
    const std::vector<Particle> &particleList, double time) {

  std::ofstream outFile;
  if (firstTimePath) {
    outFile.open(particlePathFile);
    firstTimePath = false;
  } else {
    outFile.open(particlePathFile, std::fstream::app);
  }

  if (outFile.is_open()) {
    for (auto &part : particleList) {
      outFile << boost::format("%4d ") % part.getLocation();
    }
    outFile << boost::format("%1.4e ") % time;
    outFile << "\n";
  } else {
    std::cout << "Could not open output file for simID: " << simID << "\n";
  }
}

void OutputManager::registerState(const std::vector<Particle> &particleList,
                                  double time) {

  std::ofstream outFile;
  if (firstTimeState) {
    outFile.open(stateFile);
    firstTimeState = false;
  } else {
    outFile.open(stateFile, std::fstream::app);
  }

  if (outFile.is_open()) {
    outFile << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      if (part.isAlive()) {
        outFile << boost::format("%d %d    ") % part.getType() %
                       part.getLocation();
        if (part.getType() == PType::CT) {
          outFile << boost::format("%d %d    ") % part.getType() %
                         part.getLocationCTelec();
        }
      }
    }
    outFile << "\n";
  } else {
    std::cout << "Could not open output file for simID: " << simID << "\n";
  }
}

void OutputManager::registerNumbers(const std::vector<Particle> &particleList,
                                    const Topology &topol, double time) {

  std::ofstream outFileNR;
  if (firstTimeNumber) {
    outFileNR.open(numberFile);
    firstTimeNumber = false;
  } else {
    outFileNR.open(numberFile, std::fstream::app);
  }

  std::array<int, 12> count{0};

  if (outFileNR.is_open()) {
    outFileNR << boost::format("%12.5f") % time << " ";
    for (auto &part : particleList) {
      if (part.isAlive()) {
        count[part.getType()] += 1;
        if (topol.getMolType(part.getLocation()) == MType::ben) {
          count[part.getType() + 4] += 1;
        } else if (topol.getMolType(part.getLocation()) == MType::tcne) {
          count[part.getType() + 8] += 1;
        } else {
          std::cout << "Unknown molecule type in outputmanager!\n";
          exit(EXIT_FAILURE);
        }
      }
    }

    for (auto &nr : count) {
      outFileNR << boost::format("%5d ") % nr;
    }
    outFileNR << "\n";
  }
}
