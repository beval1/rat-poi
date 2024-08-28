package org.beval;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

import java.io.File;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

public class RecombinationAnalysisParallel {

    // Constants defining the sliding window size and mismatch threshold
    private static final int WINDOW_SIZE = 100; // Size of the sliding window for sequence comparison
    private static final double THRESHOLD = 0.1; // Threshold for considering a divergence as significant
    private static String RESULT = ""; // Stores the result message
    private static final int THREADS = 4 ; //Number of threads for parallel processing

    public static void main(String[] args) {
        // Load the DNA sequences from a FASTA file
        File fastaFile = new File("src/main/resources/medium_generated_fasta_file.fasta");

        // Read the sequences from the FASTA file
        Map<String, DNASequence> dnaSequences = readFastaFile(fastaFile);

        // Analyze the sequences for recombination events using parallel processing
        List<String> recombinationEvents = analyzeRecombinationParallel(dnaSequences);

        // Uncomment the following line to print the detected recombination events
        // recombinationEvents.forEach(System.out::println);

        // Print the result message containing the time taken for analysis
        System.out.println(RESULT);
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
     * Analyzes recombination events between two DNA sequences using parallel processing.
     *
     * @param sequences A map containing DNA sequences.
     * @return A list of strings describing detected recombination events.
     */
    private static List<String> analyzeRecombinationParallel(Map<String, DNASequence> sequences) {
        // If there are fewer than two sequences, no recombination analysis can be performed
        if (sequences == null || sequences.size() < 2) return new ArrayList<>(); // Need at least two sequences to compare

        // Create a fixed thread pool executor for parallel processing
        ExecutorService executor = Executors.newFixedThreadPool(THREADS);
        List<Future<List<String>>> futures = new ArrayList<>(); // List to store Future objects for parallel tasks

        var startTime = LocalDateTime.now();

        // Extract sequences from the map for comparison
        List<DNASequence> sequenceList = new ArrayList<>(sequences.values());
        DNASequence seq1 = sequenceList.get(0); // First sequence
        DNASequence seq2 = sequenceList.get(1); // Second sequence
        int totalLength = seq1.getLength(); // Total length of the sequences

        // Dividing tasks among available processors
        int chunkSize = (totalLength - WINDOW_SIZE + THREADS - 1) / THREADS; // Determine chunk size to ensure all sequences are covered

        // Record start time for performance measurement

        // Submit tasks to the executor for parallel processing
        for (int i = 0; i < totalLength - WINDOW_SIZE; i += chunkSize) {
            final int start = i; // Start index for the chunk
            final int end = Math.min(i + chunkSize, totalLength - WINDOW_SIZE); // End index for the chunk
            // Submit a task to analyze a chunk of the sequences
            futures.add(executor.submit(() -> analyzeChunk(seq1, seq2, start, end)));
        }

        // Shutdown the executor service
        executor.shutdown();
        try {
            // Await termination of all tasks
            executor.awaitTermination(1, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            // Restore interrupted state
            Thread.currentThread().interrupt();
        }

        // Collect results from all futures
        List<String> recombinationEvents = new ArrayList<>();
        for (Future<List<String>> future : futures) {
            try {
                // Add all events detected in each chunk to the final list
                recombinationEvents.addAll(future.get());
            } catch (Exception e) {
                // Print stack trace in case of an error
                e.printStackTrace();
            }
        }

        var endTime = LocalDateTime.now();
        Duration duration = Duration.between(startTime, endTime);
        RESULT = "Time taken: " + duration.toMillis() + " milliseconds";
        return recombinationEvents; // Return the list of detected recombination events
    }

    /**
     * Analyzes a chunk of the sequences for potential recombination events.
     *
     * @param seq1  The first DNA sequence.
     * @param seq2  The second DNA sequence.
     * @param start The starting position of the chunk.
     * @param end   The ending position of the chunk.
     * @return A list of strings describing detected recombination events in the chunk.
     */
    private static List<String> analyzeChunk(DNASequence seq1, DNASequence seq2, int start, int end) {
        List<String> events = new ArrayList<>(); // List to store recombination events

        // Iterate over the chunk of sequences using a sliding window approach
        for (int i = start; i <= end; i++) {
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
        return events; // Return the list of detected recombination events in the chunk
    }
}
