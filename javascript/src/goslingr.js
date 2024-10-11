import { embed } from "gosling.js";
import { ensureWideJsonData } from "./utils"


HTMLWidgets.widget({

  name: 'goslingr',

  type: 'output',

  factory: function(el, width, height) {

      
      
      const renderValue = (x) => {
          
        // convert any long-format json data to wide-format
        const spec = ensureWideJsonData(x.spec);
    
        console.log("el!!!!")
        console.log(el)

        console.log("x!!!!")
        console.log(x)
        console.log("spec!!!!")
        console.log(spec)

        embed(document.getElementById(el.id), spec);

     
    }

    const resize = (width, height) => {

    //   TODO
    }

    return { renderValue, resize };
  }
});
