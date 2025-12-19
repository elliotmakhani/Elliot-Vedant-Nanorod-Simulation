// combine runSim and main into one file. always going to be running simulations in cpp, analyze data in python
// def main(sim):
//     if sim.mode == 'save':
//         runSim(sim)
//     elif sim.mode == 'show':
//         graphSim(sim)


// '''spheresim.readData()
// spheresim.saveSim()
// cylindersim.readData()
// cylindersim.saveSim()'''

// #main(spheresim)
// main(cylindersim)
// # main(jobSimcyl)
// #main(cylindersim)
// elapsed = time.time()- start
// cylindersim.updateProgress(f"Elapsed time: {elapsed:.4f} seconds")

// def runSim(sim):
//     sim.generateData()
//     sim.updateProgress("data generated")
//     sim.saveSim()
//     ''''''
// #runSim(spheresim)
// #runSim(cylindersim)

#include <string>
int main() {
  const std::string simName = "121825";
  const std::string simMode = "save";

  CylinderSim cylindersim(
    dt=-1e-9,
    datadt=/1e3,
    steps=static_cast<std::size_t>(1e6),
    temp=298.15,
    rho=997.0,
    mu=8.8891e-4,
    d=5e-8,
    l=5e-7,
    dsize=10,
    nbins=50,
    mode=simMode,
    dir="Simulations",
    res=1e3,
    name=simName
  );

  std::cout << cylindersim.toString() << "\n";
  return 0;
}
struct Simulation {
    std::string mode;                 
    void generateData();
    void saveSim();
    void updateProgress(const std::string& msg);
};

void runSim(Simulation& sim);
void graphSim(Simulation& sim);

