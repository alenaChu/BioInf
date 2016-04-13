import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;
import java.util.TreeSet;

/**
 * Created by Chuikina Alena on 16.3.16.
 */
public class Neighbor {
    private final char [] symbols = {'A', 'T', 'G','C'};

    public static void main (String[] args) throws IOException {
        Scanner scanner = new Scanner(System.in);
        int k = new Integer(scanner.next());
        int d = new Integer(scanner.next());
       //  String text = scanner.next();

        FileReader file = new FileReader(args[0]);
        BufferedReader br = new BufferedReader(file);
        String text = br.readLine();
        while (true) {
            try{
                text = text.concat(br.readLine());
                System.out.print(text.length() + " ");
            } catch (Exception e){
                break;
            }
        }
        System.out.println(text.length());


        Neighbor n = new Neighbor();
        for (Object i: n.frequentWordsMismatchesReverse(text, k, d))
            System.out.print(i + " ");
    }

    /** d-neighborhood of a string.
     * @param pattern
     * @param d - allowed hamming distance
     * @return The collection of strings Neighbors(Pattern, d)
     */
    public TreeSet neighbors(String pattern, int d){
        TreeSet<String> neighborhood = new TreeSet();
        neighborhood.add(pattern);
        if (d==0) {
            neighborhood.add(pattern);
            return neighborhood;
        }

        if (pattern.length()==1) {
            for (char c : symbols)  neighborhood.add(new Character(c).toString());
            return neighborhood;
        }

        String patternEnd = pattern.substring(1);
        TreeSet<String> suffixNeighbors = neighbors(patternEnd, d);
        Approx ap = new Approx();
        for (String word : suffixNeighbors){
            if (ap.hammingDistance(word, patternEnd) < d)
                for (char c: symbols){
                    String fullword = new Character(c).toString().concat(word);
                    neighborhood.add(fullword);
                }
            else {
                String fullword = pattern.substring(0,1).concat(word);
                neighborhood.add(fullword);
            }
        }
        return neighborhood;
    }

    /**  Find the most frequent k-mers with mismatches in a string.
     * @param text
     * @param k - length of pattern
     * @param d - hamming distance
     * @return All most frequent k-mers with up to d mismatches in Text
     */
    public TreeSet <String> frequentWordsWithMismatches(String text, int k, int d){
        TreeSet<String> frequentWords = new TreeSet();
        int size = (int)Math.pow(4,k);
        int [] close = new int[size];

        //close[i] set to true if this word is close to existing word in text
        PatternTransformation pt = new PatternTransformation();
        for (int i=0; i <= text.length() - k; i++)
            for (Object word: neighbors(text.substring(i,i+k),d)) {
                int index = (int) pt.patternToNumber((String) word);
                close[index]++;
            }
        //find maximum in close
        int max = 0;
        for (int i: close)
            if (i> max) max = i;
        for (int i=0; i<size; i++)
            if (close[i] == max) frequentWords.add(pt.numberToPattern(i,k));

        return frequentWords;
    }

    /** Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
     * @param text
     * @param k
     * @param d
     * @return All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Pattern)
     * over all possible k-mers
     */
    public TreeSet <String> frequentWordsMismatchesReverse(String text, int k, int d){
        TreeSet<String> frequentWords = new TreeSet();
        int size = (int)Math.pow(4,k);
        int [] close = new int[size];
        HiddenMessage hd = new HiddenMessage();

        //close[i] set to true if this word is close to existing word in text
        PatternTransformation pt = new PatternTransformation();
        for (int i=0; i <= text.length() - k; i++) {
            for (Object word : neighbors(text.substring(i, i + k), d)) {
                int index = (int) pt.patternToNumber((String) word);
                close[index]++;
            }
            for (Object word : neighbors(hd.reverseComplement(text.substring(i, i + k)), d)) {
                int index = (int) pt.patternToNumber((String) word);
                close[index]++;
            }
        }
        //find maximum in close
        int max = 0;
        for (int i: close)
            if (i> max) max = i;
        for (int i=0; i<size; i++)
            if (close[i] == max) frequentWords.add(pt.numberToPattern(i,k));

        return frequentWords;
    }

    // neighbors with one distinct symbol
    public Iterable immediateNeighbors(String pattern){
        TreeSet<String> neighborhood = new TreeSet();
        neighborhood.add(pattern);

        char [] patternArray = pattern.toCharArray();
        for (int i=0; i< patternArray.length; i++){
            char symbol = patternArray[i];
            for (char c : symbols){
                if (c == symbol) continue;
                else {
                    patternArray[i] = c;
                    neighborhood.add(new String(patternArray));
                    patternArray[i] = symbol;
                }
            }
        }
        return neighborhood;
    }

}
