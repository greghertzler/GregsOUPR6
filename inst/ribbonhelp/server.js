import * as fs from "node:fs";
import * as http from "node:http";
import * as path from "node:path";

const PORT = 8000;
const MIME_TYPES = {
  default: "application/octet-stream",
  html: "text/html; charset=UTF-16LE",
  js: "text/javascript",
  css: "text/css",
  png: "image/png",
  jpg: "image/jpg",
  gif: "image/gif",
  ico: "image/x-icon",
  svg: "image/svg+xml",
};

const STATIC_PATH = process.cwd();

const toBool = [() => true, () => false];

const prepareFile = async (pathname) => {
  const paths = [STATIC_PATH, pathname];
  if (pathname.endsWith("/")) paths.push("index.html");
  const filePath = path.join(...paths);
  const pathTraversal = !filePath.startsWith(STATIC_PATH);
  const exists = await fs.promises.access(filePath).then(...toBool);
  const found = !pathTraversal && exists;
  const streamPath = found ? filePath : STATIC_PATH + "/404.html";
  const ext = path.extname(streamPath).substring(1).toLowerCase();
  const stream = fs.createReadStream(streamPath);
  return { found, ext, stream };
};

http
  .createServer(async (req, res) => {
    const urlArray = req.url.split("?");
    const file = await prepareFile(urlArray[0]);
    const statusCode = file.found ? 200 : 404;
    const mimeType = MIME_TYPES[file.ext] || MIME_TYPES.default;
    res.writeHead(statusCode, { "Content-Type": mimeType });
    file.stream.pipe(res);
    console.log(`${req.method} ${urlArray[0]} ${statusCode}`);
  })
  .listen(PORT);

console.log(`Serving ${STATIC_PATH}\r\nListening at http://127.0.0.1:${PORT}/\r\nTo cancel, type Ctrl+C`);
