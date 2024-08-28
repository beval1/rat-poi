package org.beval;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import java.io.File;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RecombinationAnalysisSingleThread {

    // Constants defining the sliding window size and mismatch threshold
    private static final int WINDOW_SIZE = 100; // Size of the sliding window for sequence comparison
    private static final double THRESHOLD = 0.1; // Threshold for considering a divergence as significant

    public static void main(String[] args) {
        // Load the DNA sequences from a FASTA file
        File fastaFile = new File("src/main/resources/medium_generated_fasta_file.fasta");

        // Read the sequences from the FASTA file
        Map<String, DNASequence> dnaSequences = readFastaFile(fastaFile);

        // Record start time for performance measurement
        var startTime = LocalDateTime.now();

        // Analyze the sequences for recombination events
        List<String> recombinationEvents = analyzeRecombination(dnaSequences);

        // Record end time and print the time taken for the analysis
        var endTime = LocalDateTime.now();
        Duration duration = Duration.between(startTime, endTime);
        System.out.println("Time taken: " + duration.toMillis() + " milliseconds");

        // Uncomment the following line to print the detected recombination events
        // recombinationEvents.forEach(System.out::println);
    }

    /**
     * Reads a FASTA file and returns a map of sequence names to DNA sequences.
     *
     * @param file The FASTA file to read.
     * @return A map where keys are sequence names and values are DNA sequences.
     */
    private static Map<String, DNASequence> readFastaFile(File file) {
        try {
            // Read sequences from a FASTA file using BioJava's FastaReaderHelper
            return FastaReaderHelper.readFastaDNASequence(file);
        } catch (Exception e) {
            // Print stack trace in case of an error and return null
            e.printStackTrace();
            return null;
        }
    }

    /**
     * Analyzes recombination events between two DNA sequences in a sliding window manner.
     *
     * @param sequences A map containing DNA sequences.
     * @return A list of strings describing detected recombination events.
     */
    private static List<String> analyzeRecombination(Map<String, DNASequence> sequences) {
        List<String> events = new ArrayList<>(); // List to store recombination events

        // If there are fewer than two sequences, no recombination analysis can be performed
        if (sequences == null || sequences.size() < 2) return events; // Need at least two sequences to compare

        // Extract sequences from the map for comparison
        List<DNASequence> sequenceList = new ArrayList<>(sequences.values());
        DNASequence seq1 = sequenceList.get(0); // First sequence
        DNASequence seq2 = sequenceList.get(1); // Second sequence

        // Iterate over the length of the sequences using a sliding window approach
        for (int i = 0; i < seq1.getLength() - WINDOW_SIZE; i++) {
            int mismatches = 0; // Counter for mismatches in the current window

            // Compare corresponding bases in the two sequences within the current window
            for (int j = i; j < i + WINDOW_SIZE; j++) {
                // Compare bases at the current position in the two sequences
                if (seq1.getCompoundAt(j + 1).getBase().charAt(0) != seq2.getCompoundAt(j + 1).getBase().charAt(0)) {
                    mismatches++; // Increment mismatch counter if bases are different
                }
            }

            // Calculate the mismatch rate in the current window
            double mismatchRate = (double) mismatches / WINDOW_SIZE;

            // If the mismatch rate exceeds the threshold, record a recombination event
            if (mismatchRate > THRESHOLD) {
                events.add("Potential recombination detected between positions " + i + " and " + (i + WINDOW_SIZE));
            }
        }
        return events; // Return the list of detected recombination events
    }
}
