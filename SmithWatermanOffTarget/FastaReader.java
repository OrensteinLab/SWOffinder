package SmithWatermanOffTarget;
import java.nio.file.Path;
import  java.nio.file.Files;
import java.io.*;

public class FastaReader {

    public static void splitFastaToFiles(String fastaFilePath) throws IOException {
        
        String splitFolderPathStr = createTargetFolder(fastaFilePath);
        if (splitFolderPathStr == null) {
            return ;
        }

        // Open the file for reading
        BufferedReader reader = new BufferedReader(new FileReader(fastaFilePath));
        
        String line = reader.readLine();
        String header = null;
        // we will use this to append all the short sequences after the header into one string
        StringBuilder sequence = new StringBuilder();
        
        while (line != null) {
            if (line.startsWith(">")) {
                // If the line starts with ">", it is the header line
                if (header != null) {
                    // If we have already processed a sequence, save it to a file
                    saveSequenceToFile(header, sequence.toString(), splitFolderPathStr);
                    sequence = new StringBuilder();
                }
                header = line;
            } else {
                // Otherwise, it is part of the sequence
                sequence.append(line.toUpperCase());
            }
            line = reader.readLine();
        }

        // Save the last sequence to a file
        saveSequenceToFile(header, sequence.toString(), splitFolderPathStr);
        
        // Close the file
        reader.close();
    }
    
    private static String createTargetFolder(String fastaFilePath) throws IOException {
        // If there is already a split of the fasta file, then do not split
        Path filePath = Path.of(fastaFilePath);
        String fileName = filePath.getFileName().toString();
        String splitFolderPathStr = filePath.toAbsolutePath().getParent().toString() +
            "/" + fileName.substring(0, fileName.lastIndexOf('.')) + "_split/";
        Path splitFolderPath = Path.of(splitFolderPathStr);

        if (Files.exists(splitFolderPath)) {
            System.out.println(fastaFilePath + "split is alredy exists, exits without splitting.");
            return null;
        }
        else {
            // Create the folder of the split
            new File(splitFolderPathStr).mkdirs();
        }

        return splitFolderPathStr;
    }

    private static void saveSequenceToFile(String header, String sequence, String splitFolderPathStr) throws IOException {
        // Extract the sequence ID from the header
        String sequenceID = header.substring(1).trim().replaceAll("\\s", "_");
        
        // Create a new file for the sequence
        String filename = splitFolderPathStr + sequenceID + ".txt";
        FileWriter writer = new FileWriter(filename);
        
        // Write the sequence to the file
        writer.write(sequence + "\n");
        
        // Close the file
        writer.close();
        
        System.out.println("Saved sequence " + sequenceID + " to file " + filename);
    }

    public static void main(String[] args) throws IOException {
        splitFastaToFiles("genomes/hg38_only_chrs.fa");
    }
}