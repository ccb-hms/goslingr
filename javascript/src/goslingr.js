import { embed } from "gosling.js";
import { ensureWideJsonData } from "./utils"


HTMLWidgets.widget({

  name: 'goslingr',

  type: 'output',

  factory: function(el, width, height) {

      
      
      const renderValue = (x) => {
          
        // convert any long-format json data to wide-format
        const spec = ensureWideJsonData(x.spec);
        embed(el, spec);

     
    }

    const resize = (width, height) => {

    //   TODO
    }

    return { renderValue, resize };
  }
});
