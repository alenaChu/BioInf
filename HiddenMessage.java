import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by Chuikina Alena on 13.3.16.
 * Bioinformatics I
 */
public class HiddenMessage {

    public static void main(String args[]) throws IOException{
        //     Scanner scanner = new Scanner(System.in);
        FileReader file = new FileReader(args[0]);
        BufferedReader br = new BufferedReader(file);

        //  String text = scanner.next();
        String text = br.readLine();
        String arg;

        HiddenMessage hd = new HiddenMessage();

        if (args[1].compareTo("patternCount") == 0) {
            arg = br.readLine();
            int count = hd.patternCount(text, arg);
            System.out.println(count);
        }
        if (args[1].compareTo("computeFrequency") == 0) {
            arg = br.readLine();
            int [] words = hd.computingFrequency(text, new Integer(arg).intValue());
            for (int i : words) System.out.print(i + " ");
        }
        if (args[1].compareTo("frequentWord") == 0) {
            arg = br.readLine();
            Set words = hd.findingFrequentWordsBySorting(text, new Integer(arg).intValue());
            for (Object i : words) System.out.print(i + " ");
        }
        if (args[1].compareTo("reverse") == 0) {
            System.out.println(hd.reverseComplement(text));
        }
        if (args[1].compareTo("patternMatch") == 0) {
            String pattern = text;
            String genome = text;
            if (args[2] != null)
                pattern = args[2];
            else genome = br.readLine();
            for (Object i : hd.patternMatching(pattern,genome))
                System.out.print(i.toString() + " ");
        }

        br.close();
    }

    /**Implement PatternCount - subroutine in finding hidden messages in dna.
     *@param text dna String
     *@param pattern
     *@return number of Patterns in text (brute force approach)
     */
    public int patternCount(String text, String pattern){
        int k = pattern.length();
        int count =0;
        for (int i=0; i<= text.length() - k; i++)
            if (pattern.compareTo(text.substring(i, i+k)) == 0) count++;
        return count;
    }

    /**Solve the Frequent Words Problem - slow brute-force approach
     * @param text - string
     * @param k - integer = length of the substring.
     * @return All most frequent k-mers in Text.
     */
    public TreeSet frequentWords(String text, int k){
        TreeSet<String> word = new TreeSet();
        int [] count = new int[text.length() - k];
        int maxCount = 0;
        for (int i = 0; i < text.length() - k; i++){
            String pattern = text.substring(i, i + k);
            count[i] = patternCount(text, pattern);
            if (count[i] > maxCount) maxCount = count[i];
        }
        for (int i = 0; i < text.length() - k; i++)
            if (count[i] == maxCount) word.add(text.substring(i, i + k));
        return word;
    }

    /**Solve the Frequent Words Problem - using sorting
     * @param text - string
     * @param k - integer = length of the substring.
     * @return All most frequent k-mers in Text.
     */
    public TreeSet findingFrequentWordsBySorting(String text, int k){
        TreeSet<String> word = new TreeSet();
        PatternTransformation pt = new PatternTransformation();
        int maxCount = 0;
        int [] index = new int[text.length() - k];
        int [] count = new int[text.length() - k];
        for (int i = 0; i < text.length() - k; i++) {
            String pattern = text.substring(i, i + k);
            index[i] = (int)pt.patternToNumber(pattern);
            count[i] = 1;
        }
        Arrays.sort(index);
        for (int i = 1; i < index.length; i++) {
            if (index[i] == index[i-1]) count[i] = count[i-1] + 1;
            if (count[i] > maxCount) maxCount = count[i];
        }
        for (int i = 0; i<count.length; i++)
            if (count[i] == maxCount) word.add(pt.numberToPattern(i, k));
        return word;
    }

    /** the same challenge, solved using computingFrequency subroutine
     */
    public TreeSet fasterFrequentWords(String text, int k){
        TreeSet<String> word = new TreeSet();
        PatternTransformation pt = new PatternTransformation();
        int [] count = computingFrequency(text,k);
        int maxCount = 0;
        for (int i : count)
            if (i > maxCount) maxCount = i;
        for (int i = 0; i<count.length; i++)
            if (count[i] == maxCount) word.add(pt.numberToPattern(i,k));
        return word;
    }

    /** generate a frequency array of a k-mer in genome
     * @param text - genome-string
     * @param k
     * @return
     */
    private int[] computingFrequency(String text, int k){
        int pow = (int) Math.pow(4,k);
        int [] frequencyArray = new int[pow];
        PatternTransformation pt = new PatternTransformation();
        for (int i=0; i<text.length() - k+1; i++){
            String pattern = text.substring(i, i + k);
            frequencyArray[(int)pt.patternToNumber(pattern)] +=1;
        }
        return frequencyArray;
    }

    /**Reverse Complement Problem: Find the reverse complement of a DNA string.
     * @param word -A DNA string Pattern.
     * @return the reverse complement of Pattern.
     */
    public String reverseComplement(String word){
        StringBuilder temp = new StringBuilder(word);
        StringBuilder wordRev = new StringBuilder();
        for (int i = temp.length() -1 ; i>=0; i--){
            char c = temp.charAt(i);
            if (c=='T') wordRev.append('A');
            else if (c=='A') wordRev.append('T');
            else if (c=='G') wordRev.append('C');
            else if (c=='C') wordRev.append('G');
        }
        return new String(wordRev);
    }

    /** Pattern Matching Problem
     * @param pattern
     * @param text - genome-string
     * @return A collection of space-separated integers specifying all starting positions where Pattern appears
     * as a substring of Genome
     */
    public Iterable patternMatching(String pattern, String text){
        Stack <Integer> index= new Stack<>();
        int k = pattern.length();
        for (int i=0; i<= text.length() - k; i++)
            if (pattern.compareTo(text.substring(i, i+k)) == 0) index.add(new Integer(i));
        return index;
    }

    /**
     * Find patterns forming clumps in a string.
     * @param genome
     * @param k - length of pattern (k-mer)
     * @param t - minimal necessary number of k-mers in the cluster
     * @param L - length of the window
     * @return All distinct k-mers forming (L, t)-clumps in Genome
     */
    public TreeSet clumpFinding(String genome, int k, int t, int L){
        TreeSet<String> frequentPatterns = new TreeSet();
        PatternTransformation pt = new PatternTransformation();
        int size = (int) Math.pow(4, k);
        int [] clump = new int[size];

        // slide a window of length L down Genome and compute frequence array of each k-mer
        for (int i=0; i< genome.length() - L; i++){
            String text = genome.substring(i,i+L);
            int [] frequencyArray = computingFrequency(text, k);
            for (int j = 0; j< size; j++)
                if (frequencyArray[j] > t) clump[j] = 1;
        }

        for (int i=0; i< size; i++)
            if (clump[i] == 1) frequentPatterns.add(pt.numberToPattern(i,k));

        return frequentPatterns;
    }

    /**Efficient version of the previous clumpfinding function
     * Find patterns forming clumps in a string.
     * @param genome
     * @param k - length of pattern (k-mer)
     * @param t - minimal necessary number of k-mers in the cluster
     * @param L - length of the window
     * @return All distinct k-mers forming (L, t)-clumps in Genome
     */
    public TreeSet betterClumpFinding(String genome, int k, int t, int L){
        TreeSet<String> frequentPatterns = new TreeSet();
        PatternTransformation pt = new PatternTransformation();
        int size = (int) Math.pow(4, k);
        int [] clump = new int[size];

        String text = genome.substring(0,L);
        int [] frequencyArray = computingFrequency(text, k);
        for (int j = 0; j< size; j++)
            if (frequencyArray[j] > t) clump[j] = 1;
        for (int i=1; i< genome.length() - L; i++){
            String firstPattern = genome.substring(i-1,i-1+k);
            long index = pt.patternToNumber(firstPattern);
            frequencyArray[(int)index]--;
            String lastPattern = genome.substring(i+L-k,i+L);
            index = pt.patternToNumber(lastPattern);
            frequencyArray[(int)index]++;
        }

        for (int i=0; i< size; i++)
            if (clump[i] == 1) frequentPatterns.add(pt.numberToPattern(i,k));

        return frequentPatterns;
    }
}