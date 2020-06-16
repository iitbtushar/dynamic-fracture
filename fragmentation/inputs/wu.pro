log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 0";
};
  
input =
{
  file = "mesh.data";
};

// define the tolerance and max number of iteration
// for the AM solver (implemented in the ForkModule)
// AM_max_iteration = 1 => just one iteration as in Fenics code of Jian

amExp = 
{
  AM_tolerance     = 1e-4;
  AM_max_iteration = 1;
  convergence      = "Damage";
};


chain1 = 
{
  // note that: for this problem, this is just an empty file
  // as we do not need Dirichlet BCs for the phase field
  // you cannot remove the input so we need just an empty file
  input =
  {
    file = "phi.data";
  };

  model = "MP" // for parallel calculations
  {
    model = "Matrix"
    {
      matrix.type = "FEM";
      matrix.symmetric = true;
    
      model =  "VariationalDamPhase"
      {
          elements = "bulk";
          elemType = "Lagrange";
        
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };
    
          geom  = "Wu";
          thickness = 1.0;

          material  = ["mat"];

          mat = 
          {
            elements = "bulk";
            degradation=  "RandomWuDegradation"
            {
               softening="linear";
               young = 210e3;
	       ft    = 1000.;               
               Gc    = 20.;
               l0    = 1.;
               
               ft0    = 1000.;
               ftm    = 0.;
               m      = 2.;
               
               Gf0    = 20.;
               Gfm    = 2.;
            };
          };
      };
    };
  };
  
  solver = "Nonlin"
  {
      precision  = 1.0e-6;
      lineSearch = true;
      bounds     = ["b1"];

      solver = "CG"
      {
          precision = 1.0e-7;
          precon.type="ILUd";
          useThreads=true;
      };
      
      b1 = 
      {
         dofType    = "pf";
         lowerBound = 0.;
         upperBound = 1.;
      };
   };
  
  extraModules =
  {
    //modules =["energies", "view"];
    // use the following if you do not want to see the real time visualiation
    modules =["energies"];
    
    energies = "Sample"
    {
      file       = "$(CASE_NAME).dat";
      sampleWhen = "i % 2"; // every 4 steps
      header     = "Time | surface energy | strain energy | kinetic energy";
      dataSets   = [ "t","surfaceEnergy","strainEnergy", "kineticEnergy"];
    };
    

    view = "FemView"
    {
      window =
      {
        height = 300;
        width  = 900;
      };
    
      snapFile = "$(CASE_NAME)%2i.png";
      configFile  = "$(CASE_NAME).view";
    
      dataSets = [ "disp" ];
    
      disp =
      {
        type   = "Vector";
        vector = "state";
      };
    
      mesh =
      {
        //showFaces = false;
    
        plugins = "colors";
        colors =
        {
          type = "MeshColorView";
          data = "disp";
          palette   = "custom";
          autoScale = true;
        };
      };
    };
  };
};

chain2 = 
{
  input =
  {
    file = "u.data";
  };

  model = "MP"
  {
    model = "Matrix"
    {
      matrix.type = "FEM";
      matrix.symmetric = true;
    
      model       = "Multi"
      {
        models = [ "bulk", "force" ];
    
        bulk = "VariationalDamMech"
        {
          elements = "bulk";
          elemType = "Lagrange";
    
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };
        
          density   = 7850e-12;
	  E     = 210e3; //to compute CFL
    
          material = ["mat1"]; 
          mat1 =  "IsoPhaseField"
          {
            elements = "bulk";
            dim    = 2;
            state  = "PLANE_STRESS";
            thickness = 1.0;
    
            E     = 210e3;
            nu    = 0.30;
            
            degradation=  "RandomWuDegradation"
            {
               softening="linear";
               young = 210e3;
	       ft    = 1000.;               
               Gc    = 20.;
               l0    = 1.;
               
               ft0    = 1000.;
               ftm    = 0.;
               m      = 2.;
               
               Gf0    = 20.;
               Gfm    = 2.;
            };
            
            drivingForce = 
            {
              type = "Rankine";
            };
          };
        };
        
        force = "LoadScale"
        {
          // time in seconds
          scaleFunc = "return
               exp(-time/1e-4)";

          model = "LineLoad"
          {
            elements = "force";
            load     = 400.;
            shape=
            {
              type = "BLine2";
            };
          };
        };

      };
    };
  }; // end of MP model
  
  solver=  
  {
    deltaTime = 1e-8; 
  };

  extraModules = 
  {
    modules = ["vtk"];

    vtk = "VTKOutput"
    {
        prefix   = "vtu/$(CASE_NAME)";
        format   = "binary"; // binary and appended 
        elements = "bulk";
        interval = 500;
        data     = "stress | damage";
    };
  };
};



