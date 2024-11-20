package model;

import augmentedTree.Interval;

import java.util.Objects;

public class FeatureRecord implements Interval, Comparable<FeatureRecord> {

        // Fields
        private String proteinId;
        private String tag;
        private String ccdsId;
        private int start;
        private int stop;
        private double score;
        private char strand;
        private int frame;
        private String exonNumber;

        public FeatureRecord(){

        }

        // 4. start
        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        // 5. end
        public int getStop() {
            return stop;
        }

        public void setStop(int end) {
            this.stop = end;
        }

        // 6. score
        public double getScore() {
            return score;
        }

        public void setScore(double score) {
            this.score = score;
        }

        // 7. strand
        public char getStrand() {
            return strand;
        }

        public void setStrand(char strand) {
            if (strand != '-' && strand != '+' && strand != '.') {
                throw new IllegalArgumentException("Strand must be '+', '-', or '.'");
            }
            this.strand = strand;
        }

        // 8. frame
        public int getFrame() {
            return frame;
        }

        public void setFrame(int frame) {
            this.frame = frame;
        }

        // 9. protein_id
        public String getProteinId() {
            return proteinId;
        }

        public void setProteinId(String proteinId) {
            this.proteinId = proteinId;
        }


        // 12. exon_number

        public String getExonNumber() {
            return exonNumber;
        }

        public void setExonNumber(String exonNumber) {
            this.exonNumber = exonNumber;
        }

        // 14. ccds_id

        public String getCcdsId() {
            return ccdsId;
        }
        public void setCcdsId(String ccdsId) {
            this.ccdsId = ccdsId;
        }


        // 16. tag
        public String getTag() {
            return tag;
        }

        public void setTag(String tag) {
            this.tag = tag;
        }


    @Override
    public int compareTo(FeatureRecord other) {
        // First, compare by start
        int startComparison = Integer.compare(this.start, other.start);
        if (startComparison != 0) {
            return startComparison;
        }
        // If starts are equal, compare by stop
        return Integer.compare(this.stop, other.stop);
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, stop);
    }
}


