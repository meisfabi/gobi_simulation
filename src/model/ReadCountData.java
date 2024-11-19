package model;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class ReadCountData {
    private volatile ConcurrentMap<String, ConcurrentMap<String, ReadCountEntry>> readCountData;

    public Map<String, ConcurrentMap<String, ReadCountEntry>> getReadCountData() {
        if(readCountData == null) {
            synchronized(this) {
                if(readCountData == null) {
                    readCountData = new ConcurrentHashMap<>();
                }
            }
        }
        return readCountData;
    }

    public ConcurrentMap<String, ReadCountEntry> getTranscriptsForGene(String gene){
        return readCountData.get(gene);
    }
}
