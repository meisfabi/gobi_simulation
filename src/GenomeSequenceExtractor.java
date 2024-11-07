import jdk.jshell.spi.ExecutionControl;
import model.FidxEntry;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.RandomAccessFile;
import java.util.List;

public class GenomeSequenceExtractor {
    private RandomAccessFile raf;
    public GenomeSequenceExtractor(File fasta, List<FidxEntry> fidxData) throws FileNotFoundException {
        raf = new RandomAccessFile(fasta, "r");
    }

    public String getSequence(String chr, int start, int end) throws Exception{
        throw new Exception("Not Implemented");
    }


}
