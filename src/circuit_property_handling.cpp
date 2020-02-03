#include "mapper.hpp"

/**
 * initializes the circuit properties
 * 	- locations
 *  - qubits
 * if special optimization is used -> also 
 *  - depths
 *  - fidelities
 */
circuit_properties create_circuit_properties() {
    circuit_properties p;
    p.locations  = new int[nqubits];
	p.qubits     = new int[positions];
#if SPECIAL_OPT
	p.depths     = new int[positions]();
	p.fidelities = new int[positions](); 
#endif

	//Initially, no physical qubit is occupied
	for (int i = 0; i < positions; i++) {
		p.qubits[i] = -1;
	}

	//Initially, no logical qubit is mapped to a physical one
	for(unsigned i = 0; i < nqubits; i++) {
		p.locations[i] = -1;
	}

    return p;
}

int count = 0;
void adapt_circuit_properties(circuit_properties& p, const node& n) {
	delete_circuit_properties(p);
	p.locations  = n.locations;
	p.qubits     = n.qubits;
#if SPECIAL_OPT
	p.depths     = n.depths;
	p.fidelities = n.fidelities;
#endif
}

/**
 * adapts the properties of the current qubits by considering all gates of the specified layer
 */
void update_properties(circuit_properties& p, const int layer) {
#if SPECIAL_OPT
	for(std::vector<QASMparser::gate>::iterator it = layers[layer].begin(); it != layers[layer].end(); it++) {
	    QASMparser::gate g = *it;
		int pt = p.locations[g.target];
		if(g.control != -1) {
			int pc        = p.locations[g.control];
			int max_depth = std::max(p.depths[pc], p.depths[pt]) + DEPTH_GATE;
            p.depths[pc]  = max_depth;
			p.depths[pt]  = max_depth;
            p.fidelities[pt] += FIDELITY_CNOT;
            p.fidelities[pc] += FIDELITY_CNOT;

			edge e;
			e.v1 = pc;
			e.v2 = pt;
			if (graph.find(e) == graph.end()) {
				p.depths[pt]     += DEPTH_GATE    << 1; 
				p.depths[pc]     += DEPTH_GATE    << 1;
				p.fidelities[pt] += FIDELITY_GATE << 1;
				p.fidelities[pc] += FIDELITY_GATE << 1; // * 2
			}
#if USE_INITIAL_MAPPING
		} else {
#else
		} else if(pt >= 0) {
#endif
			p.depths[pt]     += DEPTH_GATE;
        	p.fidelities[pt] += FIDELITY_GATE;
		}
	}
#endif
}

/**
 * clean up the circuit properties
 */
void delete_circuit_properties(circuit_properties& p) {
    delete[] p.locations;
	delete[] p.qubits;	
#if SPECIAL_OPT
	delete[] p.depths;	
	delete[] p.fidelities;
#endif
}
