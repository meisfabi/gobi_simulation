package pooling;

import model.SimulationOutputEntry;

public class SimulationOutputFactory implements ObjectFactory<SimulationOutputEntry.Builder> {
    @Override
    public SimulationOutputEntry.Builder createObject() {
        return new SimulationOutputEntry.Builder()  ;
    }
}
