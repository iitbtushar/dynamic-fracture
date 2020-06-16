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

fork = 
{
  AM_tolerance     = 1e-4;
  AM_max_iteration = 500;
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
    
          geom = "Wu";

          material  = ["mat"];

          mat = 
          {
            elements = "bulk";
            degradation= "WuDegradation"
            {
               softening="linear";
               young = 2.4e3;
               ft    = 80.0;//Li(2016)
               Gc  = 0.2667;
               l0  = 0.0375;
            };
          };
          
          writeCracks=
          {
            file     = "$(CASE_NAME)";
            interval = 1;
            x0       = 75.;
            y0       = 35.;
            //phiMax   = 0.9;
            phiMax   = 0.85;
          };
      };
    };
  };
  
  solver = "Nonlin"
  {
      precision  = 1.0e-6;
      lineSearch = true;
      bounds     = ["b1"];

      solver = 
      {
          precision = 1.0e-7;
          type="CG";
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
      file       = "$(CASE_NAME)-energies.dat";
      sampleWhen = "i % 4"; // every 4 steps
      header     = "Time | surface energy | strain energy | kineticEnergy";
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
        models = [ "bulk", "left" ];
    
        bulk = "VariationalDamMech"
        {
          elements = "bulk";
          elemType = "Lagrange";
    
          shape =
          {
            type      = "Triangle3";
            intScheme = "Gauss1";
          };
        
          density   = 1180e-12;
    
          material = ["mat1"]; 
          mat1 = 
          {
            elements = "bulk";
            type   = "IsoPhaseField";
            dim    = 2;
            state  = "PLANE_STRESS";
    
            E     = 2.4e3;
            nu    = 0.42;
	    thickness = 1.0;
            
            degradation= "WuDegradation" 
            {
               softening="linear";
               young = 2.4e3;
               ft    = 80.0; //Li(2016)
               Gc  = 0.2667;
               l0  = 0.0375;
            };

            drivingForce = 
            {
              type = "Rankine";
            };
          };
        };
        
        left = 
        {
         type = "LoadScale";


         scaleFunc = "return
            if ( time < 0.5e-6 )
              (time*time/0.5e-6)*10.0e3
            else
             10.0e3*time
            endif";

         model =
         {
           type     = "Constraints";
           conTable = "disp";
         };
        };
      };
    };
  }; // end of MP model
  
  
  solver= "Newmark" 
  {
    deltaTime = 0.1e-7;
    solver = "Nonlin"
    {
      precision  = 1.0e-6;
      solver =
      {
          precision = 1.0e-7;
          type="CG";
          precon.type="ILUd";
          useThreads=true;
      };
    };
  };

  extraModules = 
  {
    modules = ["vtk"];

    vtk = "VTKOutput"
    {
        prefix   = "$(CASE_NAME)";
        format   = "binary"; // binary and appended 
        elements = "bulk";
        interval = 1000;
        data     = "stress|damage";
    };
  };
};



