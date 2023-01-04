/**
 * 
 */
package it.unicam.cs.asdl1819.project1;



import java.util.HashSet;
import java.util.Set;

/**
 * Implementazione dell'algoritmo di Nussinov per trovare, data una sequenza di
 * nucleotidi, una struttura secondaria senza pseudonodi che ha un numero
 * massimo di legami deboli.
 * 
 * @author Luca Tesei
 *
 */
public class NussinovFolder implements FoldingAlgorithm {

    private final String primarySequence;	// SEQUENZA DI NUCLEOTIDI DELL'RNA
    
    private int[][] matrice;			// MATRICE DI NUSSINOV.
    
    private Set<WeakBond> bondArray;	// ALL'INTERNO VERRANNO CONTENUTI TUTTI I NODI DELLA STRUTTURA OTTIMA.
    
    private boolean folded;				// VARIABILE CHE SERVIRA' PER VEDERE SE HO ESEGUITO IL METODO FOLD.
    
    private final int len ;				// LUNGHEZZA DELLA SEQUENZA PRIMARIA.
    

    /**
     * Costruisce un solver che utilizza l'algoritmo di Nussinov.
     * 
     * @param primarySequence
     *                            la sequenza di nucleotidi di cui fare il
     *                            folding
     * 
     * @throws IllegalArgumentException
     *                                      se la primarySequence contiene dei
     *                                      codici di nucleotidi sconosciuti
     * @throws NullPointerException
     *                                      se la sequenza di nucleotidi è nulla
     */
    public NussinovFolder(String primarySequence) {
        if (primarySequence == null)
            throw new NullPointerException(
                    "Tentativo di costruire un solutore Nussinov a partire da una sequenza nulla");
        String seq = primarySequence.toUpperCase().trim();
        // check bases in the primary structure - IUPAC nucleotide codes
        for (int i = 0; i < seq.length(); i++)
            switch (seq.charAt(i)) {
            case 'A':
            case 'U':
            case 'C':
            case 'G':
                break;
            default:
                throw new IllegalArgumentException(
                        "INPUT ERROR: primary structure contains an unkwnown nucleotide code at position "
                                + (i + 1));
            }

        this.primarySequence = seq;
        
        this.len = this.primarySequence.length();
        
        this.matrice = new int[len + 1][len + 1];
        
        this.bondArray = new HashSet<WeakBond>();
        
        this.folded = false;

        
        // INIZIALIZZO LA MATRICE
        for (int i = 0; i < len + 1 ; i++)
            for (int j = 0; j < len + 1; j++)
            {
            	this.matrice[i][j] = -1;
            }
    }

    public String getName() {
        return "NussinovFolder";
    }

    @Override
    public String getSequence() {
        return this.primarySequence;
    }

    @Override
    public SecondaryStructure getOneOptimalStructure() {
    	// CONTROLLO SE E' STATO ESEGUITO IL FOLD
        if(!this.folded)
        {
        	// IL FOLD NON E' STATO ESEGUITO, QUINDI LANCIO L'ECCEZIONE
        	throw new IllegalStateException("ERRORE: 	il folding sulla sequenza non è ancora stato eseguito. "); 
        }
        // ESEGUO IL TRACEBACK PER POTERMI RICAVARE UNA STRUTTURA OTTIMA ( le varie coppie che verranno trovate dal traceback verranno inserite all'interno di bounArray )
        this.traceBack(1,this.len);
        
        // VERIFICO SE IL NUMERO DEI NODI OTTENUTI DAL TRACEBACK E' EQUIVALENTE A QUELLO OTTENUTO DALL'ALGORITMO DI NUSSINOV(IN POSIZIONE [1][N]).
        if(this.matrice[1][this.len] != this.bondArray.size())	
        	throw new IllegalArgumentException("ERRORE : il folding sulla sequenza è stato eseguito in maniera errata.");
    	
    	SecondaryStructure passStructure = new SecondaryStructure(this.primarySequence,this.bondArray);
    	// VERIFICO CHE NON VI SIANO PSEUDONODI
    	if(!passStructure.isPseudoknotted())
    	{
    		// SENZA PSEUDONODI    		
    	}
    	else
    	{
    		// CI SONO SPEUDONODI
    		passStructure = null;
    		throw new IllegalArgumentException("ERRORE : Struttura secondaria con pseudonodi.");
    	}
    	
    	
        return passStructure;
    }

    @Override
    public void fold() {
    	
    	// NON ESISTE NESSUNA SEQUENZA
    	if(this.primarySequence == null)
    	{
    		throw new NullPointerException("ERRORE : non si puo eseguire il fold senza alcun dato inserito.");
    	}
    	
    	
    	//CASO BASE
    	for(int i = 1; i <= len; i++)
    	{
    		this.matrice[i][i] = 0;
    		this.matrice[i][i-1] = 0;
    	}
    	
		//CASO RICORSIVO
		
		// CICLO PRINCIPALE CHE CONTROLLA DA H = 1 ---> LEN
		for (int h = 1; h < len; h++) 
		{	
			// RIEMPIAMO LA DIAGONALE [i][j=i+h]
			for (int i = 1; i <= len - h; i++) 
			{
				// J = COLONNA 
				// I = RIGHA 
				// H = DIAGONALE
				int j = i + h;	
				for (int k = i; k < j; k++) 
				{ 
					int v=-1;
					// VERIFICO CHE I NUCLEOTIDI IN POSIZIONE (K - 1) E (J - 1) POSSANO FORMARE UNA COPPIA AMMISSIBILE
					if(this.isPaired(this.primarySequence.charAt(k - 1) , this.primarySequence.charAt(j - 1)) == 1)	// E' k - 1 e j - 1 perche all'interno dell'array l'ultimo elemento si trova
					{																							// in posizione n mentre all'interno della sequenza si trova ad n - 1
						// CALCOLO IL MASSIMO
						v = this.matrice[i][k - 1] + this.matrice[k + 1][j - 1] + 1;
						// VERIFICO SE IL MASSIMO TROVATO E' MAGGIORE DEL MASSIMO CHE SI TROVA NELLA CASELLA [I][J - 1]
						if (v > this.matrice[i][j - 1])
						{
							// INSERISCO IL MASSIMO NELLA POSIZIONE [I][J]
							this.matrice[i][j] = v;
						}
					}
					else
					{
						// I NUCLEOTIDI IN POSIZIONE (K - 1) E (J - 1) NON FORMANO UNA COPPIA AMMISIBILE, QUINDI
						// VERIFICO SE L'ELEMENTO IN POSIZIONE [I][J] E' MINORE DELL'ELEMENTO IN [I][J - 1]
						if(this.matrice[i][j] < this.matrice[i][j-1])	
						{
							// INSERISCO IL MASSIMO CHE SI TROVA NELLA CASELLA [I][J - 1] DELLA CASELLA ADIACENTE IN POSIZIONE [I][J]
							this.matrice[i][j] = this.matrice[i][j-1];
						}   		
					}
				}
			}
		}
        
        this.folded = true;
    }

    /**
     * METODO UTILIZZATO SOLAMENTE NELLA FASE DI DEBUG PER VISUALIZZARE LA
     * 		  MARICE.
     * @param matrice
     */
	private void printMatrice(int [][] matrice) {
		
		for(int i = 0; i < this.len + 1; i++)
		{
			if(i == 0)
				System.out.print("    ");
			
			System.out.print(" " + i + " |");
		}
		System.out.println();
		for(int i = 0; i < this.len + 1; i++)
		{
			if(i == 0)
				System.out.print("    ");
		
			System.out.print(" -- ");
		}
		System.out.println();
    	for(int i = 1; i <= len; i++)
    	{
    		System.out.print(" " + i + " ||");
    		for(int k = 0; k <= len; k++)
    		{
    			if(matrice[i][k] == -1)
    			{
    				System.out.print(matrice[i][k] + " |"); 
    			}
    			else
    			{
    				System.out.print(" " +matrice[i][k] + " |"); 
    			}
    		}
    		System.out.println("");
    	}
    	System.out.println("MAX = " + matrice[1][this.len]);
	}

	/**
	 * 	isPaired E' UN METODO CHE SERVE PER CONTROLLARE SE DUE NODI 
	 * 			 SI POSSONO ACCOPPIARE O MENO.
	 * 
	 * @param xi è uno dei caratteri da confrontare 
	 * @param xj è uno dei caratteri da confrontare
	 * 
	 * @return 0 se i due caratteri xi e xj non si possono legare, altrimenti 1
	 */
   private int isPaired(char xi, char xj) {
	   
   			 // NEL COSTRUTTORE TUTTA L'INTERA SEQUENZA DI CARATTERI 
   			 //		VIENE PORTATA IN 'UPPERCASE' , IO COMUNQUE HO DECISO DI METTERE
	   		 //		TUTTE LE COMBINAZIONI POSSIBILI.
	   
	      if ((xi == 'G' && xj == 'C') || (xi == 'C' && xj == 'G'))
	    	  return 1;
	      else if ((xi == 'A' && xj == 'U') || (xi == 'U' && xj == 'A'))
	    	  return 1;
	      else if ((xi == 'G' && xj == 'U') || (xi == 'U' && xj == 'G'))
	    	  return 1;
	      else if ((xi == 'g' && xj == 'c') || (xi == 'c' && xj == 'g'))
	    	  return 1;
	      else if ((xi == 'a' && xj == 'u') || (xi == 'u' && xj == 'a'))
	    	  return 1;
	      else if ((xi == 'g' && xj == 'u') || (xi == 'u' && xj == 'g'))
	    	  return 1;
	      else if ((xi == 'A' && xj == 'u') || (xi == 'u' && xj == 'A'))
	    	  return 1;
	      else if ((xi == 'a' && xj == 'U') || (xi == 'U' && xj == 'a'))
	    	  return 1;
	      else if ((xi == 'G' && xj == 'u') || (xi == 'u' && xj == 'G'))
	    	  return 1;
	      else if ((xi == 'g' && xj == 'U') || (xi == 'U' && xj == 'g'))
	    	  return 1;
	      else if ((xi == 'G' && xj == 'c') || (xi == 'c' && xj == 'G'))
	    	  return 1;
	      else if ((xi == 'g' && xj == 'C') || (xi == 'C' && xj == 'g'))
	    	  return 1;
	      else
	    	  return 0;
	   }
	

	@Override
    public boolean isFolded() {
		
		if(this.folded)
			return true;
		
        return false;
    }
    
    


    /**
     * 
     * 	traceBack serve per trovare, attraverso la matrice calcolata precedentemente a 
     *			  partire dalle posizioni (i,j),una delle strutture ottime.
     * 
     * @param i 
     * @param j
     */
	
    private void traceBack( int i, int j){
        // VERIFICO CHE L'INDICE I NON SIA MAGGIORE O UGUALE ALL'INDICE J
    	if(i>=j)
    		return ;
    	// CONTROLLO CHE LA RIGA (i+1)(j) SIA UGUALE A ALLA RIGA SOPRA
    	if((this.matrice[i + 1][j] ==  this.matrice[i][j]))
    	{
    		// SE E' UGUALE, RIESEGUO IL TRACEBACK SULLA RIGA (i+1)
    		this.traceBack(i + 1,j);
    		return;
    	}
    	// VERIFICO CHE LA COLONNA PRECENDENTE SIA UGUALE ALLA COLONNA ATTUALE
    	if((this.matrice[i][j - 1] ==  this.matrice[i][j]))
    	{	
    		// SE E' UGUALE, RIESEGUO IL TRACEBACK SULLA COLONNA PRECEDENTE
    		this.traceBack(i,j - 1);
    		return;
    	}
    	//	CONTROLLO SE IL VALORE ALL'INTERNO DELLA  DIAGONALE SX DELLA POSIZIONE CORRENTE SOMMATO CON 1(SOLAMENTE SE LE COPPIE SONO AMMISSIBILI) E' 
    	//		UGUALE AL VALORE NELLA POSIZIONE CORRENTE
    	if(((this.matrice[i + 1][j - 1]) + this.isPaired(this.primarySequence.charAt(i - 1), this.primarySequence.charAt(j - 1))) ==  this.matrice[i][j])
    	{
    		// SE E' UGUALE, MI SALVO IL NODO E RIESEGUO IL TRACEBACK SULLA DIAGONALE SX DELLA POSIZIONE CORRENTE
			this.bondArray.add(new WeakBond(i, j));
			this.traceBack(i + 1,j - 1);
    		return;
    	}
    	
		for(int k = i + 1; k < j; k++)
		{
			if((this.matrice[i][k] + this.matrice[k + 1][j]) == this.matrice[i][j])
			{
				this.traceBack(i, k);
				this.traceBack(k + 1,j);
				return;
			}
		}
    	
        return;
    }
}