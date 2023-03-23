package innovations.evolution.likelihood;


import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;


@Description("Logs innovations on the branch leading to the MRCA of a set of taxa")
public class InnovationLogger extends TreeLikelihood implements Loggable {
	public Input<TaxonSet> taxonsetInput = new Input<>("taxonset", "set of taxa defining a clade. The MRCA node of the clade is logged", Validate.REQUIRED);
	public Input<String> valueInput = new Input<>("value", "space delimited set of labels, one for each site in the alignment. Used as site label in the log file.");
	public Input<String> fromInput = new Input<>("fromStates", "states that indicate the trait is absent");
	public Input<String> toInput = new Input<>("toStates", "states that indicate the trait is present");
	
    // array of flags to indicate which taxa are in the set
    Set<String> isInTaxaSet = new LinkedHashSet<>();
    Node MRCA;
	int [] parentSample;
    int [] fromState;
    int [] toState;
    int c;
    int d;
    
    @Override
	public void initAndValidate() {
		// ensure we do not use BEAGLE
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        System.setProperty("java.only", "true");
		super.initAndValidate();
        System.setProperty("java.only", "" + forceJava);
        
        String values = valueInput.get();
		if (values != null && values.trim().length() > 0) {
			// use values as labels
			values = values.trim().replaceAll("\\s+", "\t");
			String [] strs = values.split("\t");
			if (strs.length != dataInput.get().getSiteCount()) {
				throw new IllegalArgumentException("Number of labels (" + strs.length + ") does not match amountof data (" + dataInput.get().getSiteCount() +") " + values);
			}
		}
		
		// make the to and from states as integer arrays
		parseStates();

	
        isInTaxaSet.clear();
        List<String> taxaNames = dataInput.get().getTaxaNames();
        List<String> set = taxonsetInput.get().asStringList();
        for (final String sTaxon : set) {
            final int iTaxon = taxaNames.indexOf(sTaxon);
            if (iTaxon < 0) {
                throw new IllegalArgumentException("Cannot find taxon " + sTaxon + " in data");
            }
            if (isInTaxaSet.contains(sTaxon)) {
                throw new IllegalArgumentException("Taxon " + sTaxon + " is defined multiple times, while they should be unique");
            }
            isInTaxaSet.add(sTaxon);
        }
        
//        if (fromState && taxaNames.size() == set.size()) {
//        	throw new RuntimeException("Cannot log parent of the root; either choose a different clade, or set logParent flag to false");
//        }
//        if (!fromState && ! toState) {
//        	throw new IllegalArgumentException("At least one of logMRCA and logParent should be seleceted");
//        }
	}
	
    private void parseStates() {
        // convert character string to integer array
        String fromString = fromInput.get();
        String[] fromStateStrings = fromString.split("\\s+");
        fromState = new int[fromStateStrings.length];
        for (int i = 0; i < fromStateStrings.length; i++) {
        	fromState[i] = parseInt(fromStateStrings[i], 9);
        }
        String toString = toInput.get();
        String[] toStateStrings = toString.split("\\s+");
        toState = new int[toStateStrings.length];
        for (int i = 0; i < toStateStrings.length; i++) {
        	toState[i] = parseInt(toStateStrings[i], 9);
        }
    }
    
    int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    }
    
	@Override
	public void init(PrintStream out) {
		String values = valueInput.get();
		if (values != null && values.trim().length() > 0) {
			// use values as labels
			values = values.trim().replaceAll("\\s+", "\t");
			out.append(values);
			out.append("\t");
		} else {
			int siteCount = dataInput.get().getSiteCount();
			for (int i = 0; i < siteCount; i++) {
				out.append("site"+i+"\t");
			}
		}
	}
	
	@Override
	public void log(long nSample, PrintStream out) {
		try {
			// force fresh recalculation of likelihood at this stage
			Arrays.fill(m_branchLengths, 0);
			calculateLogP();
			
			// determine the MRCA node we are going to log
            calcMRCA(treeInput.get().getRoot(), new int[1]);
            
            // sample states
            int [] sample = sample(MRCA);

            // generate output
            for (int i = 0; i < parentSample.length; i++) {
                c = parentSample[i];
                d = sample[i];
                
                if (contains(fromState, c)) { // check if c is in array fromState
                    if (contains(toState, d)) { // check if d is in array toState
                        // element at position i in array sample has a value that occurs in array toState
                        // and element at position i in array parentSample has a value that occurs in array fromState
                    	out.append(1 + "\t");
                    }else {
                    	out.append(0 + "\t");
                    }
                }else {
                	out.append(0 + "\t");
                }
            }           
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// helper method to check if an array contains a value
	public static boolean contains(int[] arr, int val) {
	    for (int i : arr) {
	        if (i == val) {
	            return true;
	        }
	    }
	    return false;
	}
	/** traverse to the root
	 * then, sample root values, and propagate back to the MRCA
	 * along the path that goes between root and MRCA
	 * @return sample
	 */
	private int[] sample(Node node) {
        int siteCount = dataInput.get().getSiteCount();
        int stateCount = dataInput.get().getMaxStateCount();
        int [] sample = new int[siteCount];

		if (node.isRoot()) {
			if (beagle != null) {
				throw new RuntimeException("BEAGLE is not supported yet");
				// m_fRootPartials = beagle.m_fRootPartials;
			}
			
			double [] p = new double[stateCount];
			double [] freqs = substitutionModel.getFrequencies();
			for (int i = 0; i < sample.length; i++) {
				int offset = stateCount * dataInput.get().getPatternIndex(i);
				for (int j = 0; j < stateCount; j++) {
					p[j] = m_fRootPartials[offset + j] * freqs[j];
				}
				sample[i] = Randomizer.randomChoicePDF(p);
			}
			
		} else {
			parentSample = sample(node.getParent());
			
			double [] p = new double[stateCount];
			double [] partials = new double[dataInput.get().getPatternCount() * stateCount * m_siteModel.getCategoryCount()];
			
			if (m_siteModel.getCategoryCount() != 1) {
				throw new RuntimeException("Gamma rate heterogeneity or proportion invariant is not supported yet");
			}
			if (beagle != null) {
				throw new RuntimeException("BEAGLE is not supported yet");
				// beagle.beagle.getPartials(arg0, arg1, arg2);
        		// getTransitionMatrix(nodeNum, probabilities);
			} else {
				likelihoodCore.getNodeMatrix(node.getNr(), 0, probabilities);
			}

			if (node.isLeaf()) { 
				if (!m_useAmbiguities.get()) {
					// leaf node values come mainly from the states.
					// only ambiguous sites are sampled
					
					int [] states = new int [dataInput.get().getPatternCount() * m_siteModel.getCategoryCount()];
					if (beagle != null) {
						throw new RuntimeException("BEAGLE is not supported yet");
						// beagle.beagle.getPartials(arg0, arg1, arg2);
		        		// getTransitionMatrix(nodeNum, probabilities);
					} else {
						likelihoodCore.getNodeStates(node.getNr(), states);
					}
					
		            for (int j = 0; j < sample.length; j++) {
		            	int childIndex = dataInput.get().getPatternIndex(j);
		            	if (states[childIndex] >= 0 && states[childIndex] < stateCount) {
		            		// copy state, if it is not ambiguous
		            		sample[j] = states[childIndex];
		            	} else {
		            		sample[j] = -1;
		            	}
		            }
				} else {
					// useAmbiguities == true
					// sample conditioned on child partials
					likelihoodCore.getNodePartials(node.getNr(), partials);

					// sample using transition matrix and parent states
		            for (int j = 0; j < sample.length; j++) {
		                int parentIndex = parentSample[j] * stateCount;
		                int childIndex = dataInput.get().getPatternIndex(j) * stateCount;
		
		                for (int i = 0; i < stateCount; i++) {
		                    p[i] = partials[childIndex + i] * probabilities[parentIndex + i];
		                }
		
						sample[j] = Randomizer.randomChoicePDF(p);
		            }
				}
			} else {
				likelihoodCore.getNodePartials(node.getNr(), partials);

				// sample using transition matrix and parent states
	            for (int j = 0; j < sample.length; j++) {
	                int parentIndex = parentSample[j] * stateCount;
	                int childIndex = dataInput.get().getPatternIndex(j) * stateCount;
	
	                for (int i = 0; i < stateCount; i++) {
	                    p[i] = partials[childIndex + i] * probabilities[parentIndex + i];
	                }
	
					sample[j] = Randomizer.randomChoicePDF(p);
	            }
            }
		}
		return sample;
	}

	@Override
	public void close(PrintStream out) {
		// nothing to do
	}
	
	
    int calcMRCA(final Node node, final int[] nTaxonCount) {
        if (node.isLeaf()) {
            nTaxonCount[0]++;
            if (isInTaxaSet.contains(node.getID())) {
                if (isInTaxaSet.size() == 1) {
                	MRCA = node;
                    return 2;
                }
                return 1;
            } else {
                return 0;
            }
        } else {
            int taxonCount = calcMRCA(node.getLeft(), nTaxonCount);
            final int nLeftTaxa = nTaxonCount[0];
            nTaxonCount[0] = 0;
            if (node.getRight() != null) {
                taxonCount += calcMRCA(node.getRight(), nTaxonCount);
                final int nRightTaxa = nTaxonCount[0];
                nTaxonCount[0] = nLeftTaxa + nRightTaxa;
                if (taxonCount == isInTaxaSet.size()) {
                	MRCA = node;
                    return taxonCount + 1;
                }
            }
            return taxonCount;
        }
    }

	
}
