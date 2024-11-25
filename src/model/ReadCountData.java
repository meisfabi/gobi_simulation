package model;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class ReadCountData {
    private volatile ConcurrentMap<String, ConcurrentMap<String, Integer>> readCountData;

    public Map<String, ConcurrentMap<String, Integer>> getReadCountData() {
        if(readCountData == null) {
            synchronized(this) {
                if(readCountData == null) {
                    readCountData = new ConcurrentHashMap<>();
                }
            }
        }
        return readCountData;
    }

    public ConcurrentMap<String, Integer> getTranscriptsForGene(String gene){
        return readCountData.get(gene);
    }
}
