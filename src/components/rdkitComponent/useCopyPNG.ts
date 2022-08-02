export async function useCopyPNG(svgText:any) {
	const convertPngBlob = async (blobSvg:any) => {
		return new Promise((resolve) => {
		 var img =new Image();
		 img.crossOrigin = "Anonymous"
		 img.src = URL.createObjectURL(blobSvg)
		 img.onload = async ()=>{
			 //@ts-ignore
			 const canvas = new OffscreenCanvas(600, 600);
			 let ctx = canvas.getContext('2d');
			 ctx.drawImage(img, 0, 0, 600, 600);
			 resolve(canvas.convertToBlob({type:'image/png', quality:1}))
		 }
	 });
	} 
	const clipboardData = {};
	const blobText = new Blob([svgText], {type:'text/plain'})
	const blobSvg = new Blob([svgText], {type:'image/svg+xml'})
	await convertPngBlob(blobSvg).then((pngBlob:any)=>{
	 clipboardData[pngBlob.type] = pngBlob
	})
	clipboardData[blobText.type] = blobText
	try {
	 console.log(clipboardData)
	 await navigator.clipboard.write([
		 new ClipboardItem(clipboardData),
	 ]);
	} catch (err:any) {
	 console.error(err.name, err.message);
	} 
}