import initRDKitModule from "@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js";
const initRDKit = (() => {
  let rdkitLoadingPromise: Promise<void>;

  return () => {
    /**
     * Utility function ensuring there's only one call made to load RDKit
     * It returns a promise with the resolved RDKit API as value on success,
     * and a rejected promise with the error on failure.
     *
     * The RDKit API is also attached to the global object on successful load.
     */
    if (!rdkitLoadingPromise) {
      rdkitLoadingPromise = new Promise((resolve, reject) => {
        initRDKitModule()
          .then((RDKit: any) => {
            window.RDKit = RDKit;
            resolve(RDKit);
          })
          .catch((e: any) => {
            reject(e);
          });
      });
    }

    return rdkitLoadingPromise;
  };
})();

export default initRDKit();
