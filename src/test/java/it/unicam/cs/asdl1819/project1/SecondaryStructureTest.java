package it.unicam.cs.asdl1819.project1;

import static org.junit.Assert.*;

import java.util.HashSet;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

public class SecondaryStructureTest {

    @Before
    public void setUp() throws Exception {
    }
    

    @Test
    public final void testHashCode() {
    	
    	SecondaryStructure first = new SecondaryStructure("AAUUGCUA");
    	
    	SecondaryStructure second = null;
    	
    	SecondaryStructure third = new SecondaryStructure("AA");
    	
    	second = first;
    	
    	assertEquals(first.hashCode(), second.hashCode());
    	
    	assertNotEquals(first.hashCode(), third.hashCode());
    	
    }

    @Test
    public final void testSecondaryStructureString() {
    	
    	SecondaryStructure structure = new SecondaryStructure("auau");
    
    	assertNotNull(structure);

    }

    @Test
    public final void testSecondaryStructureStringSetOfWeakBond() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(2, 5));
    	
    	bondArray.add(new WeakBond(3, 4));
    	
    	bondArray.add(new WeakBond(1, 6));
    	
    	SecondaryStructure structure = new SecondaryStructure("agauuu",bondArray);
    	
    	assertNotNull(structure);

    }

    @Test
    public final void testGetPrimarySequence() {
    	
    	SecondaryStructure structure = new SecondaryStructure("ccuua");
    	
    	assertEquals("CCUUA", structure.getPrimarySequence());
    	
    }

    @Test
    public final void testGetBonds() {
    	
    	int count = 0;
    	
    	NussinovFolder folder = new NussinovFolder("agu");
    	
    	folder.fold();
    	
    	SecondaryStructure structure = folder.getOneOptimalStructure();
    	
    	if(structure.getBonds().isEmpty())
    		fail("ERRORE : getBonds non funziona correttamente");
    	
    	count = structure.getBonds().size();
    	if(count != 1)
    		fail("ERRORE : getBonds non funziona correttamente");
    	
    	
    }

    @Test
    public final void testIsPseudoknotted() {
    	
    	NussinovFolder folder = new NussinovFolder("agugaugccguagcgau");
    	
    	folder.fold();
    	
    	SecondaryStructure optStructure = folder.getOneOptimalStructure();
    	
    	if(optStructure.isPseudoknotted())
    		fail("ERRORE : qualcosa ?? andato storto, sono presenti degli speudonodi.");

    }
    // BOND EMPTY
    @Test
    public final void testIsPseudoknotted1() {
    	
    	SecondaryStructure structure = new SecondaryStructure("auau");
    	
    	assertEquals(true,structure.isPseudoknotted());

    }

    // ALMENO UN LEGAME CHE HA UN ESTREMO UGUALE AD UN ALTRO LEGAME DEBOLE
    @Test
    public final void testIsPseudoknotted2() {
    	
    	SecondaryStructure structure = new SecondaryStructure("auauau");
    	
    	structure.addBond(new WeakBond(1, 4));
    	
    	structure.addBond(new WeakBond(2, 3));
    	
    	structure.addBond(new WeakBond(4, 5));
    	
    	assertEquals(true,structure.isPseudoknotted());

    }
    
    // ALMENO DUE LEGAMI SI INCROCIANO
    @Test
    public final void testIsPseudoknotted3() {
    	
    	SecondaryStructure structure = new SecondaryStructure("auauau");
    	
    	structure.addBond(new WeakBond(1, 4));
    	
    	structure.addBond(new WeakBond(3, 6));
    	
    	assertEquals(true,structure.isPseudoknotted());

    }
    
    @Test
    public final void testAddBond() {
    	
    	SecondaryStructure structure = new SecondaryStructure("auau");
    	
    	boolean status = false;
    	
    	structure.addBond(new WeakBond(2, 3));
    	
    	structure.addBond(new WeakBond(1, 4));
    	
    	status = structure.getBonds().contains(new WeakBond(2,3));
    	
    	if(!status)
    		fail("ERRORE : add bond non funziona correttamente.");
    	
    	status = structure.getBonds().contains(new WeakBond(1,4));

    	if(!status)
    		fail("ERRORE : add bond non funziona correttamente.");
    	
    	status = structure.addBond(new WeakBond(2,3));
    	
    	if(status)
    		fail("ERRORE : add bond non funziona correttamente.");
    	
    }

    @Test
    public final void testGetCardinality() {
    	
    	String sequence = "auaua";
    	
    	NussinovFolder folder = new NussinovFolder(sequence);
    	
    	folder.fold();
    	
    	SecondaryStructure structure = folder.getOneOptimalStructure();
    	
    	assertEquals(2,structure.getCardinality());
    	
    }

    @Test
    public final void testGetDotBracketNotation() {
    	
    	NussinovFolder folder = new NussinovFolder("auaua");
    	
    	folder.fold();
    	
    	SecondaryStructure structure = folder.getOneOptimalStructure();
    	
    	int status = -1;
    	
    	status = structure.getDotBracketNotation().compareTo("AUAUA\n" + 
    			" .(())");
    	
    	if(status != 0)
    		fail("ERRORE : BracketNotation non funziona correttamente");
    	
    	
    	NussinovFolder folder1 = new NussinovFolder("auauaaauugucugaucgaucgugggccccgcgcgcgcuuagacgacauucggcacauccagacaguugcagacag");
    	
    	folder1.fold();
    	
    	SecondaryStructure structure1 = folder1.getOneOptimalStructure();
    	
    	int status1 = -1;
    	
    	status1 = structure1.getDotBracketNotation().compareTo("AUAUAAAUUGUCUGAUCGAUCGUGGGCCCCGCGCGCGCUUAGACGACAUU\n"+
    			"CGGCACAUCCAGACAGUUGCAGACAG\n" + 
    			".()()(())(((((()((())((((()))(()()()()(()(.)(.)(((\n" + 
    			"()()).))..))))(())))))))..\n");

    	if(status1 != 0)
    		fail("ERRORE : BracketNotation non funziona correttamente");
    	
    }

    @SuppressWarnings("unlikely-arg-type")
	@Test
    public final void testEqualsObject() {
    	
    	SecondaryStructure firstElement = new SecondaryStructure("AUAUAU");
    	
    	SecondaryStructure secondElement = new SecondaryStructure("AUAUAU");
    	
    	NussinovFolder folder = new NussinovFolder("CCUG");
    	
    	folder.fold();
    	
    	SecondaryStructure thirdElement = folder.getOneOptimalStructure();
    	
    	if(!secondElement.equals(firstElement))
    		fail("ERRORE : equals non funziona correttamente");
    	
    	assertEquals(false,firstElement.equals(thirdElement));
    	
    	if(!firstElement.equals(firstElement))
    		fail("ERRORE : equals non funziona correttamente");
    	
    	if(!firstElement.equals(secondElement))
    		fail("ERRORE : equals non funziona correttamente");
    	
    	SecondaryStructure fourElement = null;
    	
    	assertEquals(false,firstElement.equals(fourElement));
    	
    	assertEquals(false, firstElement.equals(folder));
    	
    	firstElement.addBond(new WeakBond(1, 2));
    	
    	firstElement.addBond(new WeakBond(3, 4));
    	
    	secondElement.addBond(new WeakBond(1, 2));
    	
    	secondElement.addBond(new WeakBond(3, 4));  
    	
    	if(!secondElement.equals(firstElement))
    		fail("ERRORE : equals non funziona correttamente");
    	
    	thirdElement.addBond(new WeakBond(1, 4));
    	
    	assertEquals(false, secondElement.equals(thirdElement));
    	
    	SecondaryStructure fiveElement = new SecondaryStructure("AA");
    	
    	assertEquals(false,fiveElement.equals(firstElement));
    	
    	firstElement = new SecondaryStructure("AUAUAU");
    	
    	secondElement = new SecondaryStructure("AUAUAU");
    	
    	SecondaryStructure sixElement = new SecondaryStructure("AUAUAU");
    	
    	if(!((firstElement.equals(secondElement))&&(secondElement.equals(sixElement)&&(firstElement.equals(sixElement)))))
		{
    		fail("ERRORE : equals non funziona correttamente");
		}
    	
    	
    		
    }
    
    @Test
    public final void testToString() {
    	
    	NussinovFolder folder = new NussinovFolder("AUAU");
    	
    	folder.fold();
    	
    	SecondaryStructure el = folder.getOneOptimalStructure();
    	
    	assertEquals(el.toString(), "{(2, 3), (1, 4)}");
    	
    }
    // LA SEQUENZA PRIMARIA DI NUCLEOTIDI E' NULLA
    @Test(timeout = 10000,expected = NullPointerException.class)
    public final void testSecondaryStructureException1() {
    	
    	SecondaryStructure structure = new SecondaryStructure(null);
        	
    }
    // LA SEQUENZA PRIMARIA DI NUCLEOTIDI CONTIENE DEI CODICI DI NUCLEOTIDI SCONOSCIUTI
    @Test(timeout = 10000,expected = IllegalArgumentException.class)
    public final void testSecondaryStructureException2() {
    	
    	SecondaryStructure structure = new SecondaryStructure("AUAUAUZZTATGG");
        	
    }
    
    // LA SEQUENZA PRIMARIA DI NUCLEOTIDI E' NULLA
    @Test(timeout = 10000,expected = NullPointerException.class)
    public final void testSecondaryStructureSetOfWeakBondException1() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(2, 5));
    	
    	bondArray.add(new WeakBond(3, 4));
    	
    	bondArray.add(new WeakBond(1, 6));
    	
    	SecondaryStructure structure = new SecondaryStructure(null,bondArray);
        	
    }
    // LA SEQUENZA PRIMARIA DI NUCLEOTIDI CONTIENE DEI CODICI DI NUCLEOTIDI SCONOSCIUTI
    @Test(timeout = 10000,expected = IllegalArgumentException.class)
    public final void testSecondaryStructureSetOfWeakBondException2() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(2, 5));
    	
    	bondArray.add(new WeakBond(3, 4));
    	
    	bondArray.add(new WeakBond(1, 6));
    	
    	SecondaryStructure structure = new SecondaryStructure("AUAUAUZZTATGG",bondArray);
        	
    }
    // L'INSIEME DEI LEGAMI DEBOLI E' NULLO
    @Test(timeout = 10000,expected = NullPointerException.class)
    public final void testSecondaryStructureSetOfWeakBondException3() {
    	
    	Set<WeakBond> bondArray = null;
    	
    	SecondaryStructure structure = new SecondaryStructure("AAUUAAAUU",bondArray);
        	
    }
    // UNO DEI DUE INDICI DI UNO DEI LEGAMI DEBOLI PASSATI ESCE FUORI DAI LIMITI DELLA SEQUENZA PRIMARIA DI QUESTA STRUTTURA
    @Test(timeout = 10000,expected = IndexOutOfBoundsException.class)
    public final void testSecondaryStructureSetOfWeakBondException4() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(2, 5));
    	
    	bondArray.add(new WeakBond(3, 4));
    	
    	bondArray.add(new WeakBond(1, 6));
    	
    	SecondaryStructure structure = new SecondaryStructure("AUAU",bondArray);
        	
    }
    // ALMENO UNO DEI LEGAMI DEBOLI PASSATI CONNETTE DUE NUCLEOTIDI A FORMARE UNA COPPIA NON CONSENTITA
    @Test(timeout = 10000,expected = IllegalArgumentException.class)
    public final void testSecondaryStructureSetOfWeakBondException5() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(1, 2));
    	
    	bondArray.add(new WeakBond(3, 4));
    	
    	bondArray.add(new WeakBond(5, 6));
    	
    	SecondaryStructure structure = new SecondaryStructure("GGAUAU",bondArray);
        	
    }
    // ALMENO UNO DEI LEGAMI DEBOLI HA UN ESTREMO UGUALE AD UN ALTRO LEGAME DEBOLE.
    @Test(timeout = 10000,expected = IllegalArgumentException.class)
    public final void testSecondaryStructureSetOfWeakBondException6() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(2, 5));
    	
    	bondArray.add(new WeakBond(3, 4));
    	
    	bondArray.add(new WeakBond(1, 6));
    	
    	bondArray.add(new WeakBond(3, 5));
    	
    	SecondaryStructure structure = new SecondaryStructure("agauuu",bondArray);
        	
    }
    
    // ALMENO UNO DEI LEGAMI DEBOLI SI INCROCIANO
    @Test(timeout = 10000,expected = IllegalArgumentException.class)
    public final void testSecondaryStructureSetOfWeakBondException7() {
    	
    	Set<WeakBond> bondArray = new HashSet<WeakBond>();
    	
    	bondArray.add(new WeakBond(2, 6));
    	
    	bondArray.add(new WeakBond(3, 5));
    	
    	bondArray.add(new WeakBond(1, 4));
    	
    	SecondaryStructure structure = new SecondaryStructure("agauuu",bondArray);
        	
    }
    // IL LEGAME PASSATO E' NULLO
    @Test(timeout = 10000,expected = NullPointerException.class)
    public final void testAddBoundException() {
    	
    	SecondaryStructure structure = new SecondaryStructure("agauuu");
    	
    	structure.addBond(null);
        	
    }
    // SE INSERISCO INDICI <= 0 (DECISO ALL'INTERNO DEL COSTRUTTORE DI WEAKBOND )
    @Test(timeout = 10000,expected = IllegalArgumentException.class)	
    public final void testAddBoundException2() {
    	
    	SecondaryStructure structure = new SecondaryStructure("agauuu");
    	
    	structure.addBond(new WeakBond(0, 0));
        
    	structure.addBond(new WeakBond(-1, -1));
    	
    }
    // SE VADO FUORI LA SEQUENZA
    @Test(timeout = 10000,expected = IndexOutOfBoundsException.class)	
    public final void testAddBoundException3() {
    	
    	SecondaryStructure structure = new SecondaryStructure("agauuu");
    	
    	structure.addBond(new WeakBond(7, 9));
        
    	
    }
    // SE SI CONNETTONO DUE NUCLEOTIDI A FORMARE UNA COPPIA NON CONSENTITA
    @Test(timeout = 10000,expected = IllegalArgumentException.class)	
    public final void testAddBoundException4() {
    	
    	SecondaryStructure structure = new SecondaryStructure("agaCCu");
    	
    	structure.addBond(new WeakBond(1, 4));
    	
    }
    
    // STRUTTURA SECONDARIA CON PSEUDONODI
    @Test(timeout = 10000,expected = IllegalStateException.class)	
    public final void testGetDotBracketNotationException() {
    	
    	SecondaryStructure structure = new SecondaryStructure("auaaagua");
    	
    	structure.addBond(new WeakBond(2, 4));
    	
    	structure.addBond(new WeakBond(3, 7));
    	
    	structure.getDotBracketNotation();
    	
    }
    
    
    
}
