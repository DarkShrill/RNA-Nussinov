package it.unicam.cs.asdl1819.project1;

import static org.junit.Assert.*;



import org.junit.Before;
import org.junit.Test;

public class NussinovFolderTest {

    @Before
    public void setUp() throws Exception {
    }


    @Test
    public final void testNussinovFolder() {

        String sequence = "aUgUUgcaA";

        NussinovFolder a = new NussinovFolder(sequence);
        
        assertNotNull(a);
        
        assertEquals("AUGUUGCAA", a.getSequence());
    }

    @Test
    public final void testGetName() {
    	
        NussinovFolder a = null;
        
        a = new NussinovFolder("auau");
        
        assertEquals(a.getName(), "NussinovFolder");
    	
    }

    @Test
    public final void testGetSequence() {
    	
        NussinovFolder a = null;
        
        a = new NussinovFolder("GCUUAGCU");
        
        a.fold();
    	
        assertEquals(a.getSequence(), "GCUUAGCU");
    }

    @Test
    public final void testGetOneOptimalStructure() {
    	
        NussinovFolder a = null;
        
        a = new NussinovFolder("AAUUAUGUC");
        
        a.fold();
        
        assertNotEquals(null,a.getOneOptimalStructure());
        
    	
    }
   
    @Test
    public final void testIsFolded() {
    	
        NussinovFolder a = null;
        
        a = new NussinovFolder("AAUUAUGUC");
        
        a.fold();
        
        assertEquals(true, a.isFolded());
        	
    }
    
    public final void testIsFolded2() {
    	
        NussinovFolder a = null;
        
        a = new NussinovFolder("AAUUAUGUC");
        
        assertEquals(false, a.isFolded());
        	
    }
    
    public final void testFold() {
    	
        NussinovFolder a = null;
        
        boolean status = false;
        
        a = new NussinovFolder("AAUUAUGUC");
        
        status = a.isFolded();
        
        a.fold();
        
        status = a.isFolded();
        
        assertEquals(true, status);
        	
    }
    
    // NON ESISTE NESSUNA SEQUENZA
    @Test(timeout = 10000,expected = NullPointerException.class)
    public final void testFoldException() {
    	
        NussinovFolder a = null;
        
        a.fold();
    	
    }
    
    // TENTATIVO DI COSTRUIRE UN SOLUTORE NUSSINOV A PARTIRE DA UNA SEQUENZA NULLA
    @Test(timeout = 10000,expected = NullPointerException.class)
    public final void testNussinovFolderException1() {
    	
        NussinovFolder a = new NussinovFolder(null);
        	
    }
    // LA SEQUENZA PRIMARIA CONTINE UN NUCLEOTIDE SCONOSCIUTO 
    @Test(timeout = 10000,expected = IllegalArgumentException.class)
    public final void testNussinovFolderException2() {
    	
        NussinovFolder a = new NussinovFolder("AUAUAUTTAUAU");
        	
    }
    // IL FOLDING NON E' STATO ESEGUITO
    @Test(timeout = 10000,expected = IllegalStateException.class)
    public final void testGetOptimalStructureException() {
    	
        NussinovFolder a = new NussinovFolder("AUAAUAU");
        
        a.getOneOptimalStructure();
        	
    }
    

}
