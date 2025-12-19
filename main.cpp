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

struct Simulation {
    std::string mode;                 
    void generateData();
    void saveSim();
    void updateProgress(const std::string& msg);
};

void runSim(Simulation& sim);
void graphSim(Simulation& sim);

void main_sim(Simulation& sim) {
    if (sim.mode == "save") {
        runSim(sim);
    } else if (sim.mode == "show") {
        graphSim(sim);
    }
}


void runSim(Simulation& sim) {
    sim.generateData();
    sim.updateProgress("data generated");
    sim.saveSim();
}
