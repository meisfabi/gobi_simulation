package model;

import java.util.ArrayList;
import java.util.List;

public class SimulationOutputEntry {
    private int readId;
    private String chromosome;
    private String geneId;
    private String transcriptId;


    private String mutatedFwSeq;
    private String mutatedRvSeq;
    private long[] transcriptFwRegionVectors;
    private long[] transcriptRvRegionVectors;
    private List<Long[]> genomeFwRegionVectors;
    private List<Long[]> genomeRvRegionVectors;
    private List<Integer> fwMutationIdx;
    private List<Integer> rvMutationIdx;

    public SimulationOutputEntry(int readId, String mutatedFwSeq, String mutatedRvSeq, String chromosome, String geneId, String transcriptId, long[] transcriptFwRegionVectors, long[] transcriptRvRegionVectors, List<Long[]> genomeFwRegionVectors, List<Long[]> genomeRvRegionVectors, ArrayList<Integer> fwMutationIdx, ArrayList<Integer> rvMutationIdx) {
        this.readId = readId;
        this.chromosome = chromosome;
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.transcriptFwRegionVectors = transcriptFwRegionVectors;
        this.transcriptRvRegionVectors = transcriptRvRegionVectors;
        this.genomeFwRegionVectors = genomeFwRegionVectors;
        this.genomeRvRegionVectors = genomeRvRegionVectors;
        this.fwMutationIdx = fwMutationIdx;
        this.rvMutationIdx = rvMutationIdx;
        this.mutatedFwSeq = mutatedFwSeq;
        this.mutatedRvSeq = mutatedRvSeq;
    }

    public int getReadId() {
        return readId;
    }

    public String getChromosome() {
        return chromosome;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public long[] getTranscriptFwRegionVectors() {
        return transcriptFwRegionVectors;
    }

    public long[] getTranscriptRvRegionVectors() {
        return transcriptRvRegionVectors;
    }

    public List<Long[]> getGenomeFwRegionVectors() {
        return genomeFwRegionVectors;
    }

    public List<Long[]> getGenomeRvRegionVectors() {
        return genomeRvRegionVectors;
    }

    public List<Integer> getFwMutationIdx() {
        return fwMutationIdx;
    }

    public List<Integer> getRvMutationIdx() {
        return rvMutationIdx;
    }

    public String getMutatedFwSeq() {
        return mutatedFwSeq;
    }

    public String getMutatedRvSeq() {
        return mutatedRvSeq;
    }


}
