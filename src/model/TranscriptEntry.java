package model;

import java.util.TreeMap;
import java.util.TreeSet;

public class TranscriptEntry {


    private final TreeSet<FeatureRecord> positions = new TreeSet<>();

    public void addRecord(FeatureRecord record) {
        positions.add(record);
    }

    public TreeSet<FeatureRecord> getPositions() {
        return positions;
    }
}
