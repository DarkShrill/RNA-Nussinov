package it.unicam.cs.asdl1819.project1;


import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
//TODO:bk
/**
 * Un oggetto di questa classe rappresenta una struttura secondaria di RNA.
 * 
 * @author Luca Tesei
 *
 */
public class SecondaryStructure {

    private final String primarySequence;	// SEQUENZA DI NUCLEOTIDI DELL'RNA
    
    private final int lenPrimarySeq;		// LUNGHEZZA DELLA MIA SEQUENZA PRIMARIA

    private Set<WeakBond> bonds;			// CONTENITORE PER I MIEI LEGAMI DEBOLI
    

    /**
     * Costruisce una struttura secondaria con un insieme vuoto di legami
     * deboli.
     * 
     * @param primarySequence
     *                            la sequenza di nucleotidi
     * 
     * @throws IllegalArgumentException
     *                                      se la primarySequence contiene dei
     *                                      codici di nucleotidi sconosciuti
     * @throws NullPointerException
     *                                      se la sequenza di nucleotidi è nulla
     */
    public SecondaryStructure(String primarySequence) {
        if (primarySequence == null)
            throw new NullPointerException(
                    "Tentativo di costruire un solutore Nussinov a partire da una sequenza nulla");
        String seq = primarySequence.toUpperCase().trim();
        // check bases in the primary structure
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
        
        this.bonds = new HashSet<WeakBond>();
        
        this.lenPrimarySeq = this.primarySequence.length();
        
    }

    /**
     * Costruisce una struttura secondaria con un insieme dato di legami deboli.
     * 
     * @param primarySequence
     *                            la sequenza di nucleotidi
     * @param bonds
     *                            l'insieme dei legami deboli presenti nella
     *                            struttura
     * 
     * @throws IllegalArgumentException
     *                                       se la primarySequence contiene dei
     *                                       codici di nucleotidi sconosciuti
     * @throws NullPointerException
     *                                       se la sequenza di nucleotidi è
     *                                       nulla
     * @throws NullPointerException
     *                                       se l'insieme dei legami è nullo
     * @throws IndexOutOfBoundsException
     *                                       se almeno uno dei due indici di uno
     *                                       dei legami deboli passati esce
     *                                       fuori dai limiti della sequenza
     *                                       primaria di questa struttura
     * @throws IllegalArgumentException
     *                                       se almeno uno dei legami deboli
     *                                       passati connette due nucleotidi a
     *                                       formare una coppia non consentita.
     *                                       
     * @throws IllegalArgumentException                                     
     *                                       almeno uno dei legami deboli ha un 
	 *                                       estremo uguale ad un altro legame debole,
	 *                                       oppure due legami deboli si incociano
     * 
     */
    public SecondaryStructure(String primarySequence, Set<WeakBond> bonds) {
        if (primarySequence == null)
            throw new NullPointerException(
                    "Tentativo di costruire un solutore Nussinov a partire da una sequenza nulla");
        String seq = primarySequence.toUpperCase().trim();
        // check bases in the primary structure
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
        if(bonds == null)
        	throw new NullPointerException("ERRORE : l'insieme dei legami è nullo");
        
        this.primarySequence = seq;
        
        this.bonds = new HashSet<WeakBond>();
        
        this.lenPrimarySeq = this.primarySequence.length();
        
        
        for (WeakBond b : bonds)
        {
            this.addBond(b);
        }
        //	VERIFICO LA PRESENZA DI PSEUDONODI
    	if(this.isPseudoknotted())
    		throw new IllegalArgumentException("ERRORE : almeno uno dei legami deboli ha un estremo uguale ad un altro legame debole, oppure due legami deboli si incociano");

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
    
    /**
     * Restituisce la sequenza di nucleotidi di questa struttura secondaria.
     * 
     * @return la sequenza di nucleotidi di questa struttura secondaria
     */
    public String getPrimarySequence() {
        return this.primarySequence;
    }

    /**
     * Restituisce l'insieme dei legami deboli di questa struttura secondaria.
     * 
     * @return l'insieme dei legami deboli di questa struttura secondaria
     */
    public Set<WeakBond> getBonds() {
        return this.bonds;
    }

    
    /**
     * Determina se questa struttura contiene pseudonodi.
     * 
     * @return true, se in questa struttura ci sono almeno due legami deboli che
     *         si incrociano o se i legami sono nulli o se almeno un legame debole
     *         ha un estremo uguale ad un altro legame debole, false altrimenti
     */
    public boolean isPseudoknotted() {
    	
    	if(this.bonds.isEmpty())
    		return true;
    	
		for (WeakBond dato : this.bonds) 
		{	//	MI SALVO IL LEGAME CORRENTE
			int iPrimo = dato.getI();
			int jPrimo = dato.getJ();
			for (WeakBond temp : this.bonds) 
			{	//	SCORRO TUTTA LA LISTA DEI LEGAMI E CONFRONTO OGNI LEGAME DELLA LISTA CON
				//	IL LEGAME CORRENTE (OVVERO QUELLO SALVATO NEL FOR PRECEDENTEMENTE)
				int i = temp.getI();
				int j = temp.getJ();
				// SE MI TROVO IN PRESENZA DEL LEGAME CORRENTE( ovvero quello del ciclo for precedente) NON LO CONSIDERO.
				if(dato.equals(temp))
				{
					// NEL CASO IN CUI I DUE OGGETTI SONO UGUALI,NON FACCIO NULLA
				}
				// VEDO SE C'E' ALMENO UN LEGAME DEBOLE CHE HA UN ESTREMO UGUALE AD UN ALTRO LEGAME DEBOLE
				else if((iPrimo == i)||(jPrimo == j)||(iPrimo == j)||(jPrimo == i))
					return true;
				//VERIFICO SE CI SONO PSEUDONODI
				else if(((iPrimo < i)&&(i < j)&&(j < jPrimo))||(((i < iPrimo)&&(iPrimo < jPrimo)&&(jPrimo < j)))||((jPrimo < i))||((j < iPrimo)))
				{
					// NON FACCIO NULLA , NON CI SONO SPEUDONODI.
				}
				else
					return true;	// HO TROVATO UN PSEUDONODO.
			}
		}
		return false;
        
    }

    /**
     * Aggiunge un legame debole a questa struttura.
     * 
     * @param b
     *              il legame debole da aggiungere
     * @return true se il legame è stato aggiunto, false se era già presente
     * 
     * @throws NullPointerException
     *                                       se il legame passato è nullo
     * @throws IndexOutOfBoundsException
     *                                       se almeno uno uno dei due indici
     *                                       del legame debole passato esce
     *                                       fuori dai limiti della sequenza
     *                                       primaria di questa struttura
     * @throws IllegalArgumentException
     *                                       se il legame debole passato
     *                                       connette due nucleotidi a formare
     *                                       una coppia non consentita.
     */
    public boolean addBond(WeakBond b) {
    	
    	if(b == null)
    		throw new NullPointerException("ERRORE : il legame passato è nullo.");
    	//	VERIFICO SE ALMENO UNO DEGLI INDICI PASSATI ESCE FUORI DALLA SEQUENZA
    	if((b.getJ() <= 0 )||(b.getI() <= 0 )||(b.getJ() > this.lenPrimarySeq )||(b.getI() > this.lenPrimarySeq ))
    	{
    		throw new IndexOutOfBoundsException("ERRORE : uno dei due indici di uno dei legami deboli passati esce fuori dai limiti della sequenza primaria di questa struttura");
		} 
    	//	ALMENO UN LEGAME PASSATO CONNETTE DUE NUCLEOTIDI A FORMARE UNA COPPIA NON CONSENTITA
    	if(this.isPaired(this.primarySequence.charAt(b.getI() - 1),this.primarySequence.charAt(b.getJ() - 1)) != 1)
    	{
    		throw new IllegalArgumentException("ERRORE : almeno uno dei legami deboli passati connette due nucleotidi a formare una coppia non consentita.");
    	}
    	
    	// VERIFICO SE GIA' E' PRESENTE L'ELEMENTO B CHE VOGLIO AGGIUNGERE
    	if(!this.bonds.contains(b))
    	{
    		this.bonds.add(b);
    		
    		return true;
    	}
        return false;
    }

    /**
     * Restituisce il numero di legami deboli presenti in questa struttura.
     * 
     * @return il numero di legami deboli presenti in questa struttura
     */
    public int getCardinality() {

        return this.bonds.size();
    }



    /**
     * Restituisce una stringa contenente la rappresentazione nella notazione
     * dot-bracket di questa struttura secondaria.
     * 
     * @return una stringa contenente la rappresentazione nella notazione
     *         dot-bracket di questa struttura secondaria
     * 
     * @throws IllegalStateException
     *                                   se questa struttura secondaria contiene
     *                                   pseudonodi
     */
    public String getDotBracketNotation() {
    	
    	char[] tmpArray = new char [this.primarySequence.length() + 1];	// FACCIO PIU 1 PERCHE NON VOGLIO CONSIDERARE LA POSIZIONE 0 DELL'ARRAY. PARTO DA 1.
    	String str;
    	String strSequence;
		String string = new String();
		String stringSequence = new String();
    	
    	if(this.isPseudoknotted())
    	{
    		//LANCIO L'ECCEZIONE POICHE HO DEI PSEUDONODI
    		throw new IllegalStateException("ERRORE : questa struttura secondaria contiene pseudonodi.");
    	}
    	else
    	{
    		// SCORRO TUTTI I LEGAMI DEBOLI CHE HO E INSERISCO PER OGNI I -> '(' E PER OGNI J -> ')' ALL'INTERNO DI UN ARRAY TEMPORANEO
			for (WeakBond bond : this.bonds) {
				tmpArray[bond.getI()] = '(';
				tmpArray[bond.getJ()] = ')';
			}
			// ORA INSERISCO '.' PER TUTTE LE RESTANTI POSIZIONI NELLE QUALI NON CI SONO ELEMENTI DEL TIPO '(' O ')'
			for(char c = 1; c < this.primarySequence.length() + 1;c++)
			{
				if((tmpArray[c] != '(')&&(tmpArray[c] != ')'))
					tmpArray[c] = '.';
			}
    	}
    	str = new String(tmpArray);
    	strSequence = new String(this.primarySequence);
    	//VERIFICO SE MI TROVO IN PRESENZA DI UNA STRINGA PIU LUNGA DI 50.
    	if(str.length() >= 50)
    	{
    		int index = 0;
    		while (index < str.length()) 
    		{	// MI RICAVO UNA SCRINGA CON LUNGHEZZA MASSIMA = 50.
    		    string += (str.substring(index + 1, Math.min(index + 1 + 50 ,str.length())));	//HO MESSO PIU UNO POICHE SICCOME COME DETTO IN PRECEDENZA NON CONSIDERO
    		    string += ("\n");																// LA POSIZIONE 0 PER LA STRINGA
    		    
    		    stringSequence += (strSequence.substring(index, Math.min(index + 50 ,strSequence.length())));
    		    stringSequence += ("\n");
    		    
    		    index += 50;
    		}
    	}
    	else
    	{
    		return strSequence + "\n" + str;
    	}
        return stringSequence + string;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((bonds == null) ? 0 : bonds.hashCode());
        result = prime * result
                + ((primarySequence == null) ? 0 : primarySequence.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (!(obj instanceof SecondaryStructure))
            return false;
        SecondaryStructure other = (SecondaryStructure) obj;
        if (bonds == null) {
            if (other.bonds != null)
                return false;
        } else if (!bonds.equals(other.bonds))
            return false;
        if (primarySequence == null) {
            if (other.primarySequence != null)
                return false;
        } else if (!primarySequence.equals(other.primarySequence))
            return false;
        return true;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("{");
        Iterator<WeakBond> i = this.bonds.iterator();
        if (i.hasNext()) {
            WeakBond current = i.next();
            if(this.bonds.size() == 1)
        		sb.append(current.toString());
            while (i.hasNext()) {
                WeakBond next = i.next();
                sb.append(current.toString() + ", ");
                if (!i.hasNext()) {
                    sb.append(next.toString());
                    break;
                }
                current = next;
            }
        }

		sb.append("}");
        return sb.toString();
    }

}
