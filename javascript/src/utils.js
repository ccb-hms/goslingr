const ensureWideJsonData = (obj) => {
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        if (key === "data" && obj[key].type === "json") {
          const values = obj[key].values;
          
          // Check if values is not an array (long format)
          if (!Array.isArray(values)) {
            // Transform it to D3 format
            obj[key].values = HTMLWidgets.dataframeToD3(values);
          }
        }
  
        // If the property is an object, recursively search within it
        if (typeof obj[key] === "object" && obj[key] !== null) {
          ensureWideJsonData(obj[key]);
        }
      }
    }

    return(obj);
  };

export { ensureWideJsonData };