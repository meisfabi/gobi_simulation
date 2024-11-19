package model;

import pooling.ObjectFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

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

    public SimulationOutputEntry(int readId, String mutatedFwSeq, String mutatedRvSeq, String chromosome, String geneId, String transcriptId, long[] transcriptFwRegionVectors, long[] transcriptRvRegionVectors, List<Long[]> genomeFwRegionVectors, List<Long[]> genomeRvRegionVectors, List<Integer> fwMutationIdx, List<Integer> rvMutationIdx) {
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

    public SimulationOutputEntry(Builder builder){
        this.readId = builder.readId;
        this.chromosome = builder.chromosome;
        this.geneId = builder.geneId;
        this.transcriptId = builder.transcriptId;
        this.transcriptFwRegionVectors = builder.transcriptFwRegionVectors;
        this.transcriptRvRegionVectors = builder.transcriptRvRegionVectors;
        this.genomeFwRegionVectors = builder.genomeFwRegionVectors;
        this.genomeRvRegionVectors = builder.genomeRvRegionVectors;
        this.fwMutationIdx = builder.fwMutationIdx;
        this.rvMutationIdx = builder.rvMutationIdx;
        this.mutatedFwSeq = builder.mutatedFwSeq;
        this.mutatedRvSeq = builder.mutatedRvSeq;
    }

    public static class Builder {
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

        public Builder withReadId(int readId) {
            this.readId = readId;
            return this;
        }

        public Builder withChromosome(String chromosome) {
            this.chromosome = chromosome;
            return this;
        }

        public Builder withGeneId(String geneId) {
            this.geneId = geneId;
            return this;
        }

        public Builder withTranscriptId(String transcriptId) {
            this.transcriptId = transcriptId;
            return this;
        }

        public Builder withMutatedFwSeq(String mutatedFwSeq) {
            this.mutatedFwSeq = mutatedFwSeq;
            return this;
        }

        public Builder withMutatedRvSeq(String mutatedRvSeq) {
            this.mutatedRvSeq = mutatedRvSeq;
            return this;
        }

        public Builder withTranscriptFwRegionVectors(long[] transcriptFwRegionVectors) {
            this.transcriptFwRegionVectors = transcriptFwRegionVectors;
            return this;
        }

        public Builder withTranscriptRvRegionVectors(long[] transcriptRvRegionVectors) {
            this.transcriptRvRegionVectors = transcriptRvRegionVectors;
            return this;
        }

        public Builder withGenomeFwRegionVectors(List<Long[]> genomeFwRegionVectors){
            this.genomeFwRegionVectors = genomeFwRegionVectors;
            return this;
        }

        public Builder withFwMutationIdx(List<Integer> fwMutationIdx){
            this.fwMutationIdx = fwMutationIdx;
            return this;
        }

        public Builder withRvMutationIdx(List<Integer> rvMutationIdx){
            this.rvMutationIdx = rvMutationIdx;
            return this;
        }

        public Builder withGenomeRvRegionVectors(List<Long[]> genomeRvRegionVectors){
            this.genomeRvRegionVectors = genomeRvRegionVectors;
            return this;
        }


        public SimulationOutputEntry build() {
            return new SimulationOutputEntry(this);
        }
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
