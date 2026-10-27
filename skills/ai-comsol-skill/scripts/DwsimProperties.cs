using System;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Collections;
using System.Collections.Generic;
using System.Web.Script.Serialization;

// Direct author-maintained DWSIM thermodynamic calculator API, not a surrogate.
// This executable does not instantiate a GUI/flowsheet or save global settings.
class DwsimProperties {
  static string install;
  static Assembly Resolve(object sender, ResolveEventArgs e) {
    string path=Path.Combine(install,new AssemblyName(e.Name).Name+".dll");
    return File.Exists(path)?Assembly.LoadFrom(path):null;
  }
  [STAThread] static int Main(string[] args) {
    install=args[0]; string output=args[1];
    var report=new Dictionary<string,object>();
    var rows=new List<Dictionary<string,object>>();
    report["rows"]=rows;report["status"]="failed";
    try {
      AppDomain.CurrentDomain.AssemblyResolve+=Resolve;
      Environment.SetEnvironmentVariable("PATH",install+";"+Environment.GetEnvironmentVariable("PATH"));
      var assembly=Assembly.LoadFrom(Path.Combine(install,"DWSIM.Thermodynamics.dll"));
      dynamic calc=Activator.CreateInstance(assembly.GetType("DWSIM.Thermodynamics.CalculatorInterface.Calculator",true));
      calc.Initialize();
      calc.DisableParallelProcessing();
      string[] names=((IEnumerable)calc.GetCompoundList()).Cast<object>().Select(x=>x.ToString()).ToArray();
      string eg=names.FirstOrDefault(x=>x.Equals("Ethylene glycol",StringComparison.OrdinalIgnoreCase));
      if(eg==null)eg=names.FirstOrDefault(x=>x.ToLowerInvariant().Contains("ethylene glycol"));
      if(eg==null)throw new Exception("Ethylene glycol missing from installed compound database");
      string water=names.First(x=>x.Equals("Water",StringComparison.OrdinalIgnoreCase));
      report["dwsim_thermodynamics_version"]=assembly.GetName().Version.ToString();
      report["engine"]="DWSIM CalculatorInterface.Calculator";
      report["property_package"]="Raoult's Law";
      report["compounds"]=new[]{water,eg};
      report["constant_property_names"]=calc.GetCompoundConstPropList();
      report["temperature_dependent_property_names"]=calc.GetCompoundTDepPropList();
      report["single_phase_property_names"]=calc.GetPropList();
      if(args.Length>2 && args[2]=="pure") {
        var pure=new List<Dictionary<string,object>>();
        foreach(string compound in new[]{water,eg}) {
          var row=new Dictionary<string,object>{{"compound",compound}};
          row["molecularWeight_raw"]=calc.GetCompoundConstProp(compound,"molecularWeight");
          row["heatCapacityOfLiquid_raw_at298K"]=calc.GetCompoundTDepProp(compound,"heatCapacityOfLiquid",298.15);
          pure.Add(row);
        }
        report["pure_component_raw"]=pure;
        report["status"]="pure-property-probe";
        File.WriteAllText(Path.Combine(output,"dwsim-properties-raw.json"),new JavaScriptSerializer().Serialize(report),new System.Text.UTF8Encoding(false));
        return 0;
      }
      foreach(double temp in new[]{298.15,313.15,333.15}) {
        foreach(double xeg in new[]{0.0,0.0018,0.02,1.0}) {
          var row=new Dictionary<string,object>{{"temperature_K",temp},{"pressure_Pa_abs",101325.0},{"EG_mole_fraction",xeg}};
          foreach(string prop in new[]{"density","heatCapacityCp","viscosity","thermalConductivity"}) {
            try {
              object vals=calc.CalcProp("Raoult's Law",prop,"Mass","Liquid",new[]{water,eg},temp,101325.0,new[]{1-xeg,xeg},null,null,null,null);
              row[prop]=((IEnumerable)vals).Cast<object>().Select(x=>Convert.ToDouble(x)).ToArray();
            }catch(Exception ex){row[prop+"_error"]=ex.Message;}
          }
          rows.Add(row);
        }
      }
      report["status"]="calculated-with-property-level-status";
      report["boundary"]="neutral Water/EG reference only; excludes KOH, glycolate salts, activities, conductivity, electrochemical kinetics, membrane and selectivity";
    }catch(Exception ex){report["error"]=ex.ToString();}
    File.WriteAllText(Path.Combine(output,"dwsim-properties-raw.json"),new JavaScriptSerializer().Serialize(report),new System.Text.UTF8Encoding(false));
    return (string)report["status"]=="failed"?1:0;
  }
}
