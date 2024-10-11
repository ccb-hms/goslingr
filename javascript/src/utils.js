// wide json data is in an array (what we need):
// [{"sepal.length": 5.1, "sepal.width": 3.5}, {"sepal.length": 4.9, "sepal.width": 3}]

// dataframes are converted to long format json (smaller in size and faster to transmit over network)
// { "sepal.length": [5.1, 4.9], "sepal.width": [3.5, 3]}
// 
const ensureWideJsonData = (obj) => {
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        if (key === "data" && obj[key].type === "json") {
          const values = obj[key].values;
          
          const isWideJsonData = Array.isArray(values)
          if (!isWideJsonData) {
            // Transform it to wide format json
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