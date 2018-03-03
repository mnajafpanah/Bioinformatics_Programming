/*## CS516 (Bioinformatics Programming) Project1 
  ## Subject: Implementing the Gibbs Sampler motif finder program
  ## Author:  Mohammad Najaf-Panah
   
  
  
  #### Psudocoud ###

  #      GibbsSampler(Dna, k, t, N)
  #			randomly select k-mer Motifs = (Motif1, ..., Motif-t) in each string from Dna
  #      	BestMotifs <- Motifs
  #    	 	for j = 1 to N
  #				i = RANDOM(t)
  #				Profile <- profile matrix formed from all strings in Motifs except for Motif-i
  #				Motif-i <- Profile randomly generated k-mer in the i-th sequence
  #				if Score(Motifs) < Score(BestMorifs)
  #					BestMorifs <- Motifs
  #			return BestMorifs
  #           
*/


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.stream.Collectors;



public class GibbsSampler {
    
	private static final int LARGE_NUMBER = 100000;
	private static  int [] [] baseCounts;
	
	/**
     * A method to read a FASTA format file which contains sequences and their descriptions
     * @param file
     * @return a list of DNA sequences
     */
    	public static List<String> readFASTAFile(String file) {
    	
    		List<String> desc= new ArrayList<>();
    		List<String> seq = new ArrayList<>();
	
    		try{
    			BufferedReader input = new BufferedReader( new FileReader( file ) );
    			StringBuffer buffer = new StringBuffer();
    			String line = input.readLine();
     
    				if(line == null)
    					throw new IOException( file + "Error: is an empty file! Try another file." );
    					//System.exit(0);
     
    				if( line.charAt(0) != '>' )
    					throw new IOException( "Error: First line of " + file + " should start with '>'." );
    				
    				else
    					desc.add(line); 
        
        for( line = input.readLine().trim(); line != null; line = input.readLine() ) {
            if( line.length()>0 && line.charAt(0) == '>' ) {
            		seq.add(buffer.toString());
            		buffer = new StringBuffer();
            		desc.add(line);
            } else  
            		buffer.append(line.trim());
        	} // end for-loop   
        
        if(buffer.length() != 0)
        		seq.add(buffer.toString());
    		}catch(IOException e) {
    			System.out.println("Error when reading ... " + file);
    			e.printStackTrace();
    		}
	
	return seq;
 
    } // end readFASTAFile
	
    	
    	/**
         * A method which returns lexicographical index related to each Nucleotide Base.
         * @param baseLetter    Nucleotide Base
         * @return  lexicographical index related to the letter
         */
        public static long letterIndex(char baseLetter) {
            
        	switch (baseLetter) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                case 'a': return 0;
                case 'c': return 1;
                case 'g': return 2;
                case 't': return 3;
            }
        	
            return -1;
            
        } // end symbolIndex
        
        
        /**
         * A method to calculate the motif count  matrix using a given list of motifs.  
         * @param motifs    given motifs 
         * @return          motif count matrix which is a 4 by (motif length) matrix that gives each Base counts at each position
         */
        public static double[][] MotifCountMatrix(List<String> motifs) {
            
        		int numMotifs = motifs.size();
            
        		int motifLength = motifs.get(0).length();
            
        		double[][] count = new double[4][motifLength];
            
        		for (String motif : motifs) {
                
        			for (int j = 0; j < motifLength; j++) {
                    
        				char base = motif.charAt(j);
                    
        				int baseIndex = (int) letterIndex(base);
                    
        				count[baseIndex][j]++;
                }
            
        		} // end for-loop

            return count;
            
        } // end of MotifCountMatrix
        
     
        /**
         * A method which given a list of motifs and returns the probability profile matrix which is a 4 x (motif length) matrix.
         * @param motifs    a list of given motifs
         * @return          profile matrix which consist of probability of each base (A, C, G, T) in each motif
         */
        public static double[][] ProfileGenetator(List<String> motifs) {
            
        		int numMotifs = motifs.size();
            int motifLength = motifs.get(0).length(); // value t
            double[][] profile = new double[4][motifLength]; // it's a 4 x (motif length) matrix
            
            // go over the profile matrix
            for (String motif : motifs) {
                
            		for (int j = 0; j < motifLength; j++) {
                    
            			char base = motif.charAt(j);
                    int baseIndex = (int) letterIndex(base);
                    profile[baseIndex][j] += 1.0 / numMotifs;
                }
            } // end for-loop

            return profile;
            
        } // end ProfileGenetator
        
        
    
    /**
     * A method to calculate the probability of a k-mer given a probability profile matrix.
     * @param kmerPattern   given k-mer
     * @param profile    probability profile matrix
     * @return probability of pattern (k-mer)
     */
    public static double KmerProbability(String kmerPattern, double[][] profile) {
        
    		double probability = 1.0;
        
    		for (int i = 0; i < kmerPattern.length(); i++) {
            
    			char base = kmerPattern.charAt(i);
            
    			probability *= profile[(int) letterIndex(base)][i];
        }
    		
        return probability;
        
    } // end 
    
    /**
     * A method which given a profile matrix and a sequence. Then the method will randomly select one of its k-mers.
     * @param profile   profile matrix
     * @param text      it's a string sequence of DNA
     * @return          random k-mer generated from profile
     */
    public static String KmerMotifGenerator(double[][] profile, String sequence) {
        int k = profile[0].length;
        List<Double> probabilities = new ArrayList<>();
        
        for (int i = 0; i < sequence.length() - k; i++) {
            String pattern = sequence.substring(i, i + k);
            probabilities.add(KmerProbability(pattern, profile)); // collecting the probability of each k-mer in a given sequence
        } // end for-loop
        
        double totalProbability = probabilities.stream().mapToDouble(Double::doubleValue).sum();
        
        Random random = new Random();
        double randomProbability = random.nextDouble()*totalProbability;
        double runningSum = 0.;
        
        for (int i = 0; i < sequence.length() - k; i++) {
            runningSum += probabilities.get(i);
            if (runningSum > randomProbability) {
                return sequence.substring(i, i + k);
            }
        } // end for-loop
        
        return sequence.substring(sequence.length() - k);
    
    } // end KmerMotifGenerator
    
    
    /**
     * A method for applying Laplace's Rule of succession to form profile from Motif1 to Motif(i-1). So, a small number (called pseudocounts) 
     * such as 1 is added to each entry in the unnormalized matrix.
     * @param motifs    motifs matrix
     * @return profile matrix after applying the Laplace's Rule
     */
    public static double[][] LaplaceProfile(List<String> motifs) {
        
    		int numMotifs = motifs.size(); 
        int motifLength = motifs.get(0).length(); // value t
        double[][] profile = new double[4][motifLength];

        for (int baseIndex = 0; baseIndex <= 3; baseIndex++) {
            for (int j = 0; j < motifLength; j++) {
                profile[baseIndex][j] += 1.0 / (numMotifs + 4.0);
            }
        } // end for-loop

        for (String motif : motifs) {
            for (int j = 0; j < motifLength; j++) {
                char base = motif.charAt(j);
                int baseIndex = (int) letterIndex(base);
                profile[baseIndex][j] += 1.0 / (numMotifs + 4.0);
            }
        } // end for-loop

        return profile;
    
    } // end LaplaceProfile
    
    
    
    /**
     * A method to calculate hamming distance score among all motifs compare with consensus. Actually, the score is 
     * the total number of mismatches from the consensus string (the most common base at each position) over all motifs.
     * @param motifs given motifs 
     * @return  Hamming distance score of motifs
     */
    public static int ScoreGenerator(List<String> motifs) {
        
    		int HDscore = 0;
        int numMotifs = motifs.size();
        double[][] profile = ProfileGenetator(motifs);
        
        // go over profile matrix
        for (int j = 0; j < profile[0].length; j++) {
            
        		List<Double> BaseCounts = new ArrayList<>(Arrays.asList(profile[0][j], profile[1][j], profile[2][j], profile[3][j]));
            
            HDscore += numMotifs * (1.0 - Collections.max(BaseCounts)); // calculating the # of mismatches in each iteration
        
        } // end for-loop
        
        return HDscore;
        
    } // end ScoreGenerator
    

    /**
     * A method which given k, t, N and a list of sequences and returns the best motifs with the lowest score.
     * @param Dna  a matrix which holds the sequences
     * @param k    length of k-mers
     * @param t    number of sequences for sampling
     * @param N    number of sampling or searches
     * @return   a list of k-mer motifs with lowest score (best motifs)
     */
    public static List<String> GibbsSampler(List<String> Dna, int k, int t, int N) {
        
    		int textLength = Dna.get(0).length();
        Random random = new Random();
        
        // Get random list of motifs
        List<String> motifs = Dna.stream().map(s -> {
            int index = random.nextInt(textLength - k);
            return s.substring(index, index + k);
        }).collect(Collectors.toList());
        
        List<String> BestMotifs = new ArrayList<>(motifs);
        
        for (int j = 0; j < N; j++) {
        	
            int i = random.nextInt(t);
            List<String> sublistTexts = new ArrayList<>(motifs);
            sublistTexts.remove(i);
            double[][] profile = LaplaceProfile(sublistTexts);
            motifs.set(i, KmerMotifGenerator(profile, Dna.get(i)));
            if (ScoreGenerator(motifs) < ScoreGenerator(BestMotifs)) {
                BestMotifs = motifs;
            }
            
        } // end for-loop
        
        return BestMotifs;
        
    } // end GibbsSampler

    /**
     * Gibbs Sampler method which returns the best k-mer motifs with the lowest score after running the GibbsSampler R times (searches).
     * @param Dna    a list of sequences
     * @param k      length of each k-mers
     * @param t      number of sequences 
     * @param N      number of samplings or searches
     * @param R      number of searches (Runs)
     * @return       a list of k-mer motifs with lowest score
     */
    public static List<String> GibbsSamplerViaMultiSearches(List<String> Dna, int k, int t, int N, int R) {
       
    		int bestScore = LARGE_NUMBER;
        List<String> bestMotifs = GibbsSampler(Dna, k, t, N);
        
        for (int run = 0; run < R - 1; run++) {
            List<String> motifs = GibbsSampler(Dna, k, t, N);
            int score = ScoreGenerator(motifs);
            if (score < bestScore) {
                bestMotifs = motifs;
                bestScore = score;
            }
            
        }
        
        //System.out.println(bestScore);
        
        return bestMotifs;
        
    } // end of GibbsSamplerViaMultiSearches
    
    /**
     * A method to generate the consensus sequence with a given motif matrix
     * @param motifs
     */
    public static String consensusGenerator(List<String> motifs) {
    	
    		Double maxNum = 0.0;
    		String answer = "";
    		double [] [] countMatrix = MotifCountMatrix(motifs);
    		
    		// Go over count matrix
    		for(int j = 0; j < countMatrix[0].length; j++) {
    	
    		List<Double> BaseCounts = new ArrayList<>(Arrays.asList(countMatrix[0][j], countMatrix[1][j], countMatrix[2][j], countMatrix[3][j]));
        maxNum = Collections.max(BaseCounts);
        if(maxNum == countMatrix[0][j])
        		answer += "A";
        else if(maxNum == countMatrix[1][j])
        		answer += "C";
        else if(maxNum == countMatrix[2][j])
        		answer += "G";
        else if(maxNum == countMatrix[3][j])
        		answer += "T";
        
    		}
    		
    		return answer;
    		
    } // end of consensusGenerator

    /**
     * A method to calculate logarithm base 2.
     * @param base
     * @param num
     * @return log base 2 of num
     */
    public static double log(int base, Double num) {
        
    		return Math.log(num) / Math.log(base);
    
    }
    
    /**
     * A method to calculate the entropy of the motif profile
     * @param motifs
     * @param k
     * @return entropy of the motif profile
     */
    public static double entropyCalculator(List<String> motifs, int k) {
    		//DecimalFormat df = new DecimalFormat("#.##");
    		double [] [] profilMatrix = ProfileGenetator(motifs);
        double entropy = 0.00;
        
        // Go over profile matrix to get each probability
        for(int i = 0; i < 4; i++) {
        		for(int j = 0; j < k; j++) {
        			if(profilMatrix[i][j] != 0)
        			entropy += profilMatrix[i][j] * (log(2, profilMatrix[i][j]));
        		}
        } // end for-loop
        
        return  -(entropy);
    
    } // end entropyCalculator
    
    
    //###########################################       MAIN        ###################################/////
    //###############
    //#####
    /**
     * This is the main program which execute the test cases and also demonstrate the results for a given FASTA file.
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
    	
    		/////////////##############     TEST CASE ONE     ###############/////////////
    		String fileName1 = "";
		if (args.length>0) 
			fileName1 = args[0];
		else {
			System.out.print("TEST CASE ONE for a 4-mer\n");
			System.out.print("Please enter the 1st FASTA file name (test1.fasta): "); // ask users to enter the name of "text file"
			Scanner sc1 = new Scanner(System.in);
			fileName1 = sc1.next();
		}
		List<String> Dna1 = readFASTAFile(fileName1); // read test1.fasta file 
		System.out.print("By givening arbitrary 4-mer motifs: ");
		int 	k1 = 4; 
    		int N1 = 200;
    		int R1 = 400;
    		int t1 = Dna1.size();
    		List<String> motifs1 = GibbsSamplerViaMultiSearches(Dna1, k1, t1, N1, R1);
    		List<String> answer1 = Arrays.asList("tttc", "tttt", "ttct", "tttt");
    		
    		// Check Count Matrix
    		double [] [] countMatrix1 = MotifCountMatrix(motifs1);
    		List<Double> Motifscount = new ArrayList<>();
    		for(int i = 0; i < 4; i++) {
    			for(int j = 0; j < 4; j++) {
    			Motifscount.add(countMatrix1[j][i]);
    			}
    		}
    		List<Double> answer2 = Arrays.asList(0.0, 0.0, 0.0,4.0, 0.0, 0.0, 0.0,4.0, 0.0, 1.0, 0.0,3.0, 0.0, 1.0, 0.0,3.0);
    		
    		
    		// Check Profile Matrix
    		double [] [] profilMatrix1 = ProfileGenetator(motifs1);
    		List<Double> MotifsProfile = new ArrayList<>();
    		for(int i = 0; i < 4; i++) {
    			for(int j = 0; j < 4; j++) { 
    				MotifsProfile.add(profilMatrix1[j][i]);
    			}
    		}
    		List<Double> answer3 = Arrays.asList(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 0.75, 0.0, 0.25, 0.0, 0.75);
    		
    		// Check Consensus Sequence
    		String consensus1 = consensusGenerator(motifs1);
    		String answer4 = "TTTT";
    	
    		// Check Hamming distance score of the motif
    		int HammingDis = ScoreGenerator(motifs1);
    		int answer5 = 2;	
    		
    		// Check Entropy of the motif profile
    		DecimalFormat df = new DecimalFormat("0.00");
    		df.format(entropyCalculator(motifs1, k1));
        String answer6 = "1.62";
    		
    		if(motifs1.equals(answer1) && Motifscount.equals(answer2) && MotifsProfile.equals(answer3) && consensus1.equals(answer4) && HammingDis==answer5 && ((df.format(entropyCalculator(motifs1, k1))).equals(answer6))) {
    			System.out.print("\n\n###################################################");
    			System.out.print("\n\n########## Test ONE Successfully Passed.   ########");
    			System.out.print("\n\n###################################################");
    			System.out.print("\n\nThis is my EXPECTED motifs: " + answer1);
    			System.out.print("\nThis is my OBSERVED motifs: ");
    	        System.out.println("\n" + motifs1.stream().collect(Collectors.joining("\n")));
    			System.out.print("\nThis is my EXPECTED motif Count: " + answer2);
    			System.out.print("\nThis is my OBSERVED motif Count: " + Motifscount.toString());
    			System.out.print("\nThis is my EXPECTED motif Profile: " + answer3);
    			System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile.toString());
    			System.out.print("\nThis is my EXPECTED Consensus: " + answer4);
    			System.out.print("\nThis is my OBSERVED Consensus: " + consensus1);
    			System.out.print("\nThis is my EXPECTED Hamming distance score: " + answer5);
    			System.out.print("\nThis is my OBSERVED Hamming distance score: " + HammingDis);
    			System.out.print("\nThis is my EXPECTED Entropy of the motif profile: " + answer6);
    			System.out.print("\nThis is my OBSERVED Entropy of the motif profile: " + (df.format(entropyCalculator(motifs1, k1))));
    		
    		} else {
    			
    			System.out.print("\n\nError: Unfortunately Test 1 Failed!\n");
    			System.out.print("\nThis is my EXPECTED motifs: " + answer1);
    			System.out.print("\nThis is my OBSERVED motifs: " + motifs1);
    	        //System.out.println("\n" + motifs1.stream().collect(Collectors.joining("\n")));
    			//System.out.print("\nThis is my OBSERVED motif Count: " + countMatrix1 +"\n");
    			System.out.print("\nThis is my EXPECTED motif Count: " + answer2);
    			System.out.print("\nThis is my OBSERVED motif Count: " + Motifscount);
    			System.out.print("\nThis is my EXPECTED motif Profile: " + answer3);
    			System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile.toString());
    			System.out.print("\nThis is my EXPECTED Consensus: " + answer4);
    			System.out.print("\nThis is my OBSERVED Consensus: " + consensus1);
    			System.out.print("\nThis is my EXPECTED Hamming distance score: " + answer5);
    			System.out.print("\nThis is my OBSERVED Hamming distance score: " + HammingDis);
    			System.out.print("\nThis is my EXPECTED Entropy of the motif profile: " + answer6);
    			System.out.print("\nThis is my OBSERVED Entropy of the motif profile: " + (df.format(entropyCalculator(motifs1, k1))));
    			System.out.print("\n\n");
    		}
     
        
    		
		/////////////##############     TEST CASE TWO    ###############/////////////
    		String fileName2 = "";
		if (args.length>0) 
			fileName2 = args[0];
		else {
			System.out.print("\n\n\nTEST CASE TWO for a 4-mer");
			System.out.print("\n\n");
			System.out.println("Please enter the 2nd FASTA file name (test2.fasta): "); // ask users to enter the name of "text file"
			Scanner sc2 = new Scanner(System.in);
			fileName2 = sc2.next();
		}
		List<String> Dna2 = readFASTAFile(fileName2); // read test2.fasta file 
		System.out.print("By givening arbitrary 4-mer motifs: ");
		int 	k2 = 4; 
    		int N2 = 200;
    		int R2 = 728;
    		int t2 = Dna2.size();
    		List<String> motifs2 = GibbsSamplerViaMultiSearches(Dna2, k2, t2, N2, R2);
    		List<String> answer2_1 = Arrays.asList("gtgt", "gtat", "gcat", "ctat");
    		
    		// Check Count Matrix
    		double [] [] countMatrix2 = MotifCountMatrix(motifs2);
    		List<Double> Motifscount2 = new ArrayList<>();
    		for(int i = 0; i < 4; i++) {
    			for(int j = 0; j < 4; j++) {
    				Motifscount2.add(countMatrix2[j][i]);
    			}
    		}
    		List<Double> answer2_2 = Arrays.asList(0.0, 1.0, 3.0, 0.0, 0.0, 1.0, 0.0, 3.0, 3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 4.0);
    		
    		
    		// Check Profile Matrix
    		double [] [] profilMatrix2 = ProfileGenetator(motifs2);
    		List<Double> MotifsProfile2 = new ArrayList<>();
    		for(int i = 0; i < 4; i++) {
    			for(int j = 0; j < 4; j++) { 
    				MotifsProfile2.add(profilMatrix2[j][i]);
    			}
    		}
    		List<Double> answer2_3 = Arrays.asList(0.0, 0.25, 0.75, 0.0, 0.0, 0.25, 0.0, 0.75, 0.75, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 1.0);
    		// Check Consensus Sequence
    		String consensus2 = consensusGenerator(motifs2);
    		String answer2_4 = "GTAT";
    	
    		// Check Hamming distance score of the motif
    		int HammingDis2 = ScoreGenerator(motifs2);
    		int answer2_5 = 3;	
    		
    		// Check Entropy of the motif profile
    		DecimalFormat df2 = new DecimalFormat("0.00");
    		df2.format(entropyCalculator(motifs2, k2));
        String answer2_6 = "2.43";
    		
    		if(motifs2.equals(answer2_1) && Motifscount2.equals(answer2_2) && MotifsProfile2.equals(answer2_3) && consensus2.equals(answer2_4) && HammingDis2==answer2_5 && ((df2.format(entropyCalculator(motifs2, k2))).equals(answer2_6))) {
    			System.out.print("\n\n###################################################");
    			System.out.print("\n\n########## Test Two Successfully Passed.   ########");
    			System.out.print("\n\n###################################################");
    			System.out.print("\n\nThis is my EXPECTED motifs: " + answer2_1);
    			System.out.print("\nThis is my OBSERVED motifs: ");
    	        System.out.println("\n" + motifs2.stream().collect(Collectors.joining("\n")));
    			System.out.print("\nThis is my EXPECTED motif Count: " + answer2_2);
    			System.out.print("\nThis is my OBSERVED motif Count: " + Motifscount2.toString());
    			System.out.print("\nThis is my EXPECTED motif Profile: " + answer2_3);
    			System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile2.toString());
    			System.out.print("\nThis is my EXPECTED Consensus: " + answer2_4);
    			System.out.print("\nThis is my OBSERVED Consensus: " + consensus2);
    			System.out.print("\nThis is my EXPECTED Hamming distance score: " + answer2_5);
    			System.out.print("\nThis is my OBSERVED Hamming distance score: " + HammingDis2);
    			System.out.print("\nThis is my EXPECTED Entropy of the motif profile: " + answer2_6);
    			System.out.print("\nThis is my OBSERVED Entropy of the motif profile: " + (df2.format(entropyCalculator(motifs2, k2))));
    		
    		} else {
    			
    			System.out.print("\n\nError: Unfortunately Test 2 Failed!\n");
    			System.out.print("\nThis is my EXPECTED motifs: " + answer2_1);
    			System.out.print("\nThis is my OBSERVED motifs: " + motifs2);
    	        System.out.println("\n" + motifs2.stream().collect(Collectors.joining("\n")));
    			System.out.print("\nThis is my OBSERVED motif Count: " + countMatrix1 +"\n");
    			System.out.print("\nThis is my EXPECTED motif Count: " + answer2_2);
    			System.out.print("\nThis is my OBSERVED motif Count: " + Motifscount2);
    			System.out.print("\nThis is my EXPECTED motif Profile: " + answer2_3);
    			System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile2.toString());
    			System.out.print("\nThis is my EXPECTED Consensus: " + answer2_4);
    			System.out.print("\nThis is my OBSERVED Consensus: " + consensus2);
    			System.out.print("\nThis is my EXPECTED Hamming distance score: " + answer2_5);
    			System.out.print("\nThis is my OBSERVED Hamming distance score: " + HammingDis2);
    			System.out.print("\nThis is my EXPECTED Entropy of the motif profile: " + answer2_6);
    			System.out.print("\nThis is my OBSERVED Entropy of the motif profile: " + (df2.format(entropyCalculator(motifs2, k2))));
    			System.out.print("\n\n");
    		}
     
    		
    		/////////////##############     TEST CASE Three    ###############/////////////
    		String fileName3 = "";
		if (args.length>0) 
			fileName3 = args[0];
		else {
			System.out.print("\n\n\nTEST CASE THREE for a 5-mer");
			System.out.print("\n\n");
			System.out.println("Please enter the 3rd FASTA file name (test3.fasta): "); // ask users to enter the name of "text file"
			Scanner sc3 = new Scanner(System.in);
			fileName3 = sc3.next();
		}
		List<String> Dna3 = readFASTAFile(fileName3); // read test3.fasta file 
		System.out.print("By givening arbitrary 5-mer motifs: ");
		int 	k3 = 5; 
    		int N3 = 200;
    		int R3 = 400;
    		int t3 = Dna3.size();
    		List<String> motifs3 = GibbsSamplerViaMultiSearches(Dna3, k3, t3, N3, R3);
    		List<String> answer3_1 = Arrays.asList("cctta", "tctgt", "cgttc", "cccta", "cgtca");
    		//List<String> answer3_1 = Arrays.asList("cctt", "tctg", "cgtt", "ccct", "cgtc");
    		
    		// Check Count Matrix
    		double [] [] countMatrix3 = MotifCountMatrix(motifs3);
    		List<Double> Motifscount3 = new ArrayList<>();
    		for(int i = 0; i < 5; i++) {
    			for(int j = 0; j < 4; j++) {
    	
    				Motifscount3.add(countMatrix3[j][i]);
    			}
    		}
    		List<Double> answer3_2 = Arrays.asList(0.0, 4.0, 0.0, 1.0, 0.0, 3.0, 2.0, 0.0, 0.0, 1.0, 0.0, 4.0, 0.0, 1.0, 1.0, 3.0, 3.0, 1.0, 0.0, 1.0);
    		
    		
    		
    		// Check Profile Matrix
    		double [] [] profilMatrix3 = ProfileGenetator(motifs3);
    		List<Double> MotifsProfile3 = new ArrayList<>();
    		for(int i = 0; i < 5; i++) {
    			for(int j = 0; j < 4; j++) { 
    				MotifsProfile3.add(profilMatrix3[j][i]);
    			}
    		}
    		System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile3.toString());
    		List<Double> answer3_3 = Arrays.asList(0.0, 0.8, 0.0, 0.2, 0.0, 0.6000000000000001, 0.4, 0.0, 0.0, 0.2, 0.0, 0.8, 0.0, 0.2, 0.2, 0.6000000000000001, 0.6000000000000001, 0.2, 0.0, 0.2);
    		
    		// Check Consensus Sequence
    		String consensus3 = consensusGenerator(motifs3);
    		String answer3_4 = "CCTTA";
    	
    		// Check Hamming distance score of the motif
    		int HammingDis3 = ScoreGenerator(motifs3);
    		int answer3_5 = 3;	
    		
    		// Check Entropy of the motif profile
    		DecimalFormat df3 = new DecimalFormat("0.00");
    		df3.format(entropyCalculator(motifs3, k3));
        String answer3_6 = "5.16";
        
        //
    		if(motifs3.equals(answer3_1)  && Motifscount3.equals(answer3_2)  && MotifsProfile3.equals(answer3_3) && consensus3.equals(answer3_4) && HammingDis3==answer3_5 && ((df3.format(entropyCalculator(motifs3, k3))).equals(answer3_6))) {
    			System.out.print("\n\n###################################################");
    			System.out.print("\n\n########## Test Three Successfully Passed. ########");
    			System.out.print("\n\n###################################################");
    			System.out.print("\n\nThis is my EXPECTED motifs: " + answer3_1);
    			System.out.print("\nThis is my OBSERVED motifs: ");
    	        System.out.println("\n" + motifs3.stream().collect(Collectors.joining("\n")));
    			System.out.print("\nThis is my EXPECTED motif Count: " + answer3_2);
    			System.out.print("\nThis is my OBSERVED motif Count: " + Motifscount3.toString());
    			System.out.print("\nThis is my EXPECTED motif Profile: " + answer3_3);
    			System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile3.toString());
    			System.out.print("\nThis is my EXPECTED Consensus: " + answer3_4);
    			System.out.print("\nThis is my OBSERVED Consensus: " + consensus3);
    			System.out.print("\nThis is my EXPECTED Hamming distance score: " + answer3_5);
    			System.out.print("\nThis is my OBSERVED Hamming distance score: " + HammingDis3);
    			System.out.print("\nThis is my EXPECTED Entropy of the motif profile: " + answer3_6);
    			System.out.print("\nThis is my OBSERVED Entropy of the motif profile: " + (df3.format(entropyCalculator(motifs3, k3))));
    		
    		} else {
    			
    			System.out.print("\n\nError: Unfortunately Test 3 Failed!\n");
    			System.out.print("\nThis is my EXPECTED motifs: " + answer3_1);
    			System.out.print("\nThis is my OBSERVED motifs: " + motifs3);
    	        System.out.println("\n" + motifs1.stream().collect(Collectors.joining("\n")));
    			System.out.print("\nThis is my OBSERVED motif Count: " + countMatrix1 +"\n");
    			System.out.print("\nThis is my EXPECTED motif Count: " + answer3_2);
    			System.out.print("\nThis is my OBSERVED motif Count: " + Motifscount3);
    			System.out.print("\nThis is my EXPECTED motif Profile: " + answer3_3);
    			System.out.print("\nThis is my OBSERVED motif Profile: " + MotifsProfile3.toString());
    			System.out.print("\nThis is my EXPECTED Consensus: " + answer3_4);
    			System.out.print("\nThis is my OBSERVED Consensus: " + consensus3);
    			System.out.print("\nThis is my EXPECTED Hamming distance score: " + answer3_5);
    			System.out.print("\nThis is my OBSERVED Hamming distance score: " + HammingDis3);
    			System.out.print("\nThis is my EXPECTED Entropy of the motif profile: " + answer3_6);
    			System.out.print("\nThis is my OBSERVED Entropy of the motif profile: " + (df3.format(entropyCalculator(motifs3, k3))));
    		}
       
        
        //########################///////       START PROGRAM        //////######################
    		//############
    		//####
    		System.out.print("\n\n#################################################################");
    		System.out.print("\n\n###########   GIBBS SAMPLER FINDING MOTIF PROGRAM   #############");
		System.out.print("\n\n#################################################################");
    		long startTime = System.nanoTime();
        
    		// THIS IS THE PROGRAM TO PROCESS THE TWO BIG FILES
    		String fileName = "";
    		
    		if (args.length>0) 
    			fileName = args[0];
    		else {
    			System.out.print("\n\nPlease enter your FASTA file name (50seq1000bp.fasta OR 500seq100bp.fasta): "); // ask users to enter the name of "text file"
    			Scanner sc = new Scanner(System.in);
    			fileName = sc.next();
    		}
    		
    		List<String> Dna = readFASTAFile(fileName);
    		
    		System.out.print("Please enter the size of k (k-mer): "); // ask users to enter the value of "k"
    		Scanner scan = new Scanner(System.in);
    		int k = scan.nextInt();
    		int t = Dna.size(); 
    		int N = 200; // number of samplings or searches 
    		int R = 200; // number of run for Gibss Sampler	
    		List<String> motifs = GibbsSamplerViaMultiSearches(Dna, k, t, N, R);
    		System.out.print("\nThe motif from each sequence is: ");
        System.out.println("\n" + motifs.stream().collect(Collectors.joining("\n")));
        
        
        // Print the motif count matrix
        double [] [] countMatrix = MotifCountMatrix(motifs);
        System.out.print("\nThe motif count matrix is: " + "\n");
        for(int i = 0; i < 4; i++) {
        		for(int j = 0; j < k; j++) {
        			System.out.print((int)countMatrix[i][j]);
        			System.out.print("    ");
        		}
        System.out.print("\n");
        }
        
        // Print the motif profile matrix
        double [] [] profilMatrix = ProfileGenetator(motifs);
        System.out.println("\nThe motif profile matrix is: " + "\n");
        for(int i = 0; i < 4; i++) {
        		for(int j = 0; j < k; j++) {
        			System.out.printf("%.2f", profilMatrix[i][j]);
        			System.out.print("  ");
        		}
        System.out.print("\n");
        }
        
        
        // Print consensus sequence
        String consensus = consensusGenerator(motifs);
        System.out.print("\n\nThe consensus sequence is: " + consensus); 
   
        
        // Print hamming distance score of the motif
        System.out.print("\n\nThe Hamming distance score of the bestMotifs is: " + ScoreGenerator(motifs));
        
        
        // Print Entropy of the motif profile
        Double entropy = entropyCalculator(motifs, k);
        System.out.print("\n\nThe Entropy of the motif profile is: " + entropy); 

        long endTime   = System.nanoTime();
		long totalTime = endTime - startTime;
		System.out.println("\n\nRuninng time: " + totalTime + " nanoseconds");
        
        
    } // END OF MAIN
    
} // END OF CLASS