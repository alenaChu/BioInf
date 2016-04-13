import java.util.Random;
import java.util.Scanner;
import java.util.Stack;
import java.util.TreeSet;

/**
 * Created by Chuikina Alena on 19.3.16.
 * BioInformatics I, weeks 3 and 4
 */
public class MedianString {
    private final char [] symbols = {'A','C','G','T'};

    public static void main(String[] args) {
        // put your code here
        Scanner scanner = new Scanner(System.in);
        //   String pattern = scanner.next();
        int k = scanner.nextInt();
        int t = scanner.nextInt();
        int N = scanner.nextInt();
        String[] dna = new String[t];
        for (int i=0; i< t; i++)
            dna[i] = scanner.next();
        MedianString ms = new MedianString();

        String [] bestMotifs = new String [t];
        for (int i=0; i< 100; i++) {
            String[] motifs = ms.gibbsSampler(dna, k, t, N);
            System.out.println(ms.computeScore(motifs));
            if (bestMotifs[0] == null || ms.computeScore(bestMotifs) > ms.computeScore(motifs))
                for (int ind = 0; ind < motifs.length; ind++) bestMotifs[ind] = motifs[ind];
        }
        for (Object motif: bestMotifs)
            System.out.println(motif);
        System.out.println(ms.computeScore(bestMotifs));
        //System.out.print(ms.profileKMer(text,k,profile));


        //  System.out.print(ms.medianString(dna,k) + " ");
        //for (Object i : ms.createKMers(k))
        //    System.out.print(i + " ");
    }

    /** brute force algorithm for finding motifs in dna string-arrays
     *  Motif - pattern that occure in all given strings (with some mismatches)
     * @param dna - collection of strings
     * @param k
     * @param d - max number of mismatches
     * @return All (k, d)-motifs in Dna
     */
    public TreeSet<String> motifEnumerator(Stack<String> dna, int k, int d){
        TreeSet <String> patterns = new TreeSet<>();
        Neighbor nb = new Neighbor();
        for (String string: dna){
            for (int i=0; i<= string.length() - k; i++){
                String pattern = string.substring(i,i+k);
                TreeSet neighbors = nb.neighbors(pattern,d);
                for (Object pattern1: neighbors){
                    TreeSet<String> neighbors1 = nb.neighbors((String)pattern1,d);
                    boolean flag = true;
                    for (String string1: dna) {
                        int j;
                        for (j =0; j<= string1.length() - k; j++){
                            String word = string1.substring(j,j+k);
                            if (neighbors1.contains(word)) {
                                break;
                            }
                        }
                        if (j == string1.length()-k + 1) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) patterns.add((String)pattern1);
                }
            }
        }
        return patterns;
    }

    /**MedianString
     * @param dna - collection of strings
     * @param k
     * @return Any k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern.
     */
    public String medianString(String [] dna, int k){
        String Median = new String();
        TreeSet kMers = createKMers(k);
        int dist = Integer.MAX_VALUE;
        Approx ap = new Approx();
        for (Object kmer: kMers){
            int d = ap.d((String)kmer,dna);
            if (dist > d){
                dist = d;
                Median = (String)kmer;
            }
        }

        return Median;
    }
    //subroutine in medianString search
    private TreeSet createKMers (int k){
        TreeSet kMers = new TreeSet();
        if (k==1) {
            for (char c : symbols)  kMers.add(new Character(c).toString());
            return kMers;
        }

        TreeSet<String> suffixKMer = createKMers(k - 1);
        for (String word : suffixKMer){
            for (char c: symbols){
                String fullword = new Character(c).toString().concat(word);
                kMers.add(fullword);
            }
        }
        return kMers;
    }

    /** find a Profile-most probable k-mer in a string.
     * @param text
     * @param k
     * @param profile - 4 × k matrix with probabilities of each symbol in a row
     * @return  A Profile-most probable k-mer in Text
     */
    public String profileKMer(String text, int k, double [][] profile){
        String pattern = text.substring(0, k);
        double max = 0;
        for (int i=0; i<= text.length()-k; i++){
            String pattern1 = text.substring(i, i + k);
            double prob = probability(pattern1, profile);
            if (prob > max) {
                max = prob;
                pattern = pattern1;
            }
        }
        return pattern;
    }

    /** Motif search - greedy approach - at each step chose motif with best probability
     * @param dna - collection  of strings
     * @param k
     * @param t - length if dna
     * @return A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna,k,t).
     */
    public String [] greedyMotifSearch(String [] dna, int k, int t){
        String [] bestMotifs = new String[dna.length];
        for (int i=0;i<dna.length; i++)
            bestMotifs[i] = dna[i].substring(0,k);
        double minScore = computeScore(bestMotifs);

        for (int i=0; i<= dna[0].length() - k; i++){
            String [] motifs = new String[dna.length];
            motifs[0] = dna[0].substring(i,i+k);
            for (int j=1; j<dna.length; j++){
                double [][] profile = profileMatrix(motifs,k);
                motifs[j] = profileKMer(dna[j],k,profile);
            }
            double score = computeScore(motifs);
            if (score < minScore) {
                bestMotifs = motifs;
                minScore = score;
            }
        }
        return bestMotifs;
    }

    /** run randomized algorithm for motif search 1000 times and save best result
     * @param dna
     * @param k
     * @return  A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k) 1,000 times
     */
    public String[] randomizedMotifSearch(String[] dna, int k){
        String [] globalBestMotifs = new String[dna.length];
        double globalMinScore = Integer.MAX_VALUE;

        for (int iter=0; iter< 1000; iter++) {
            String [] bestMotifs = randomSearch1(dna, k);
            double minScore = computeScore(bestMotifs);
            if (globalMinScore > computeScore(bestMotifs) ) {
                globalBestMotifs = bestMotifs;
                globalMinScore = minScore;
            }
        }
        return globalBestMotifs;
    }

    /** another randomized algorithm for motif search  - GibbsSampler
     * @param dna - array
     * @param k
     * @param t - number of strings in dna
     * @param N = 20 - how many times to start gibbsSampler
     * @return
     */
    public String[] gibbsSampler(String[] dna, int k, int t, int N){
        String [] bestMotifs = new String[dna.length];
        Random randomGenerator = new Random();
        for (int j = 0; j < dna.length; j++){
            int rand = randomGenerator.nextInt(dna[0].length() - k +1);
            bestMotifs[j] = dna[j].substring(rand,rand+k);
        }

        String [] motifs = new String[dna.length];
        for (int i=0; i< motifs.length; i++) motifs[i] = bestMotifs[i];

        for (int j=0; j<N; j++){
            int i = (int)(t*Math.random());
            motifs[i] = "";
            double [][]profile = profileMatrix(motifs, k);
            int index = random(probabilities(dna[i], profile, k));
            motifs[i] = dna[i].substring(index,index+k);
            if (computeScore(motifs) < computeScore(bestMotifs))
                for (int ind=0; ind< motifs.length; ind++) bestMotifs[ind] = motifs[ind];
        }
        return bestMotifs;
    }
    /***************************************************************************
     * Helper functions
     ***************************************************************************/
    // subroutine for randomized motif search
    private String [] randomSearch1(String[] dna, int k){
        String [] bestMotifs = new String[dna.length];
        for (int j = 0; j < dna.length; j++){
            int rand = (int)  (  Math.random() * (dna[0].length() - k) );
            bestMotifs[j] = dna[j].substring(rand,rand+k);
        }

        double minScore = computeScore(bestMotifs);
        String[] motifs = bestMotifs;
        while(true){
            double[][] profile = profileMatrix(motifs,k);
            for (int j = 0; j< dna.length; j++)
                motifs[j] = profileKMer(dna[j], k, profile);

            double score = computeScore(motifs);
            if (score < minScore) {
                bestMotifs = motifs;
                minScore = score;
            }
            else return bestMotifs;
        }
    }

    // find probability of each k-mer in a string - used in random function
    private double[] probabilities(String text, double [][] profile, int k){
        double [] prob = new double [text.length() - k+1];
        double sum =0 ;
        for (int i=0; i< prob.length; i++) {
            prob[i] = probability(text.substring(i, i + k), profile);
            sum+=prob[i];
        }
        for (int i=0; i< prob.length; i++)
            prob[i] /= sum;
        return prob;
    }

    // random number generator for gibbsSampler
    private int random(double[] prob){
        double i = Math.random();
        double sum = 0;
        int j;
        for (j=0; j<prob.length; j++) {
            sum += prob[j];
            if (i < sum) return j;
        }
        return j;
    }

    // create 4*k matrix with probabilities of each of 4 symbols at the position
    private double [][] profileMatrix(String [] motifs, int k){
        // array can be not full
        int numWords = motifs.length;
        for (int i = numWords-1; i>=0; i--)
            if (motifs[i] == null) numWords--;

        double [][] profile = new double[4][k];
        char [][] motifMatrix = new char[numWords][k];
        for (int i=0; i<motifs.length; i++)
            if (motifs[i] != null && motifs[i].length() > 0)
                motifMatrix[i] = motifs[i].toCharArray();

        for (int i=0; i<motifs.length; i++)
            if (motifs[i]!=null && motifs[i].length() > 0)
                for (int j=0; j< k; j++)
                    profile[getRow(motifMatrix[i][j])][j]++;

        for (int i=0; i<4; i++)
            for (int j=0; j< k; j++)
                profile[i][j] = (profile[i][j]+1)/(numWords+4);
        return profile;
    }

    // score the set of motifs
    private double computeScore(String[] motifs){
        int numWords = motifs.length;
        for (int i = numWords-1; i>0; i--)
            if (motifs[i] == null) numWords--;

        int wordLength = motifs[0].length();
        double [][] profile = new double[4][wordLength];
        char [][] motifMatrix = new char[numWords][wordLength];
        for (int i=0; i<numWords; i++)
            motifMatrix[i] = motifs[i].toCharArray();

        for (int i=0; i<numWords; i++)
            for (int j=0; j< wordLength; j++)
                profile[getRow(motifMatrix[i][j])][j]++;

        double score = 0;
        for (int i=0; i < profile[0].length; i++) {
            // find maximum in column
            double max = 0;
            for (int j = 0; j < 4; j++) {
                score += profile[j][i];
                if (profile[j][i] > max)
                    max = profile[j][i];
            }
            score -= max;
        }

        return score;
    }

    //another way to score the motif
    private double computeEntropy(String[] motifs, int k){
        double [][] profile = profileMatrix(motifs, k);

        double entropy = 0;
        for (int i=0; i < profile.length; i++)
            for (int j = 0; j< profile[0].length; j++)
                if (profile[i][j] > 0)
                    entropy -= profile[i][j]*Math.log(profile[i][j])/Math.log(2);

        return entropy;
    }

    // for profile matrix
    private int getRow(char c){
        if (c == 'A' || c == 'a')return 0;
        if (c == 'C' || c == 'c') return 1;
        if (c == 'G' || c == 'g') return 2;
        if (c == 'T' || c == 't') return 3;
        else return -1;
    }

    // the probability of a word according to profile matrix
    private double probability(String word, double [][]profile){
        double prob = 1;
        char [] wordArray = word.toCharArray();
        for (int i=0; i<word.length(); i++) {
            int row = getRow(wordArray[i]);
            prob *= profile[row][i];
        }
        return prob;
    }


}
