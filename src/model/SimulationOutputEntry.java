package model;

import java.util.List;
import java.util.Set;

public class SimulationOutputEntry {

    private int readId;
    private int[] transcriptFwRegionVectors;
    private int[] transcriptRvRegionVectors;
    private List<int[]> genomeFwRegionVectors;
    private List<int[]> genomeRvRegionVectors;
    private Set<Integer> fwMutationIdx;
    private Set<Integer> rvMutationIdx;

    public SimulationOutputEntry(Builder builder) {
        this.readId = builder.readId;
        this.transcriptFwRegionVectors = builder.transcriptFwRegionVectors;
        this.transcriptRvRegionVectors = builder.transcriptRvRegionVectors;
        this.genomeFwRegionVectors = builder.genomeFwRegionVectors;
        this.genomeRvRegionVectors = builder.genomeRvRegionVectors;
        this.fwMutationIdx = builder.fwMutationIdx;
        this.rvMutationIdx = builder.rvMutationIdx;
    }

    public static class Builder {
        private int readId;
        private int[] transcriptFwRegionVectors;
        private int[] transcriptRvRegionVectors;
        private List<int[]> genomeFwRegionVectors;
        private List<int[]> genomeRvRegionVectors;
        private Set<Integer> fwMutationIdx;
        private Set<Integer> rvMutationIdx;

        public Builder withReadId(int readId) {
            this.readId = readId;
            return this;
        }

        public Builder withTranscriptFwRegionVectors(int[] transcriptFwRegionVectors) {
            this.transcriptFwRegionVectors = transcriptFwRegionVectors;
            return this;
        }

        public Builder withTranscriptRvRegionVectors(int[] transcriptRvRegionVectors) {
            this.transcriptRvRegionVectors = transcriptRvRegionVectors;
            return this;
        }

        public Builder withGenomeFwRegionVectors(List<int[]> genomeFwRegionVectors) {
            this.genomeFwRegionVectors = genomeFwRegionVectors;
            return this;
        }

        public Builder withFwMutationIdx(Set<Integer> fwMutationIdx) {
            this.fwMutationIdx = fwMutationIdx;
            return this;
        }

        public Builder withRvMutationIdx(Set<Integer> rvMutationIdx) {
            this.rvMutationIdx = rvMutationIdx;
            return this;
        }

        public Builder withGenomeRvRegionVectors(List<int[]> genomeRvRegionVectors) {
            this.genomeRvRegionVectors = genomeRvRegionVectors;
            return this;
        }


        public SimulationOutputEntry build() {
            return new SimulationOutputEntry(this);
        }
    }

    public int[] getTranscriptFwRegionVectors() {
        return transcriptFwRegionVectors;
    }

    public int[] getTranscriptRvRegionVectors() {
        return transcriptRvRegionVectors;
    }

    public List<int[]> getGenomeFwRegionVectors() {
        return genomeFwRegionVectors;
    }

    public List<int[]> getGenomeRvRegionVectors() {
        return genomeRvRegionVectors;
    }

    public Set<Integer> getFwMutationIdx() {
        return fwMutationIdx;
    }

    public Set<Integer> getRvMutationIdx() {
        return rvMutationIdx;
    }

    public int getReadId() {
        return readId;
    }
}
